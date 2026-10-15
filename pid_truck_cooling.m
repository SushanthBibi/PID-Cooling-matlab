%% PID Thermal Control for Heavy-Duty Trucks
% Keeps the engine coolant temperature T_c at 90 °C.
% Transfer functions, step response, P vs PI vs PID
% (The Simulink model is built by hand from the same numbers)
% Sushanth Bibi (24BME0501), Ganishk Mohanraj (24BME0500)

clc;
clear;
close all;

%% Physical data (from published sources)
% Volvo D13C460 fact sheet: 12.8 L heavy-duty truck engine
P_rated = 338;     % rated power [kW] (460 hp)
m_dry   = 1121;    % engine dry weight [kg]
V_cool  = 38;      % total coolant volume [L]
% Material properties (Engineering ToolBox)
cp_cool = 3.5;     % 50/50 water-glycol specific heat [kJ/kg K]
rho     = 1.04;    % 50/50 water-glycol density [kg/L]
cp_iron = 0.46;    % cast iron specific heat [kJ/kg K]
cp_al   = 0.91;    % aluminium specific heat [kJ/kg K]
% Heavy-duty cooling module data (Worsztynowicz, 2017)
Q_max   = 150;     % heat to coolant at max power [kW] (paper: 133-167)
T_max   = 105;     % coolant temperature at that point [°C]
T_hot   = 45;      % ambient temperature of that test [°C]

%% Assumptions (not from a source, stated openly)
T_ref = 90;        % set-point [°C]
T_amb = 25;        % normal ambient temperature [°C]
L0    = 30;        % normal engine load [%]
f_m   = 0.5;       % share of engine mass that heats up with the coolant
V_eng = 26;        % coolant inside the engine [L] (rest is in radiator + hoses)
m_rad = 30;        % radiator core mass [kg]
flow  = 250;       % pump flow [L/min] (inside the paper's range)
% Heat to coolant is proportional to load, and the cooling
% conductance grows linearly with the command y (0-3).

%% Derive the plant constants (linearised energy balance)
% Heat capacities [kJ/K]:  C = mass x specific heat
C_c = V_eng*rho*cp_cool + f_m*m_dry*cp_iron;        % engine + its coolant
C_r = (V_cool - V_eng)*rho*cp_cool + m_rad*cp_al;   % radiator + its coolant
mcp = flow/60*rho*cp_cool;                          % coolant flow [kW/K]

% Cooling conductance H = heat removed per K above ambient [kW/K]
H_max = Q_max/(T_max - T_hot);          % full cooling (y = 3)
q_L   = Q_max/100;                      % heat to coolant per % load [kW/%]
H0    = q_L*L0/(T_ref - T_amb);         % needed at 90 °C and 30 % load
y0    = 3*H0/H_max;                     % cooling command at that point

% Time constants [s] and gains
tau_c  = C_c/H0;                        % engine:   C / H
tau_r  = C_r/mcp;                       % radiator: C / (m*cp)
K_cool = (H_max/3)*(T_ref - T_amb)/H0;  % fall in T_c per unit of y [K]
K_load = q_L/H0;                        % rise in T_c per % load [K/%]

disp('--- Derived plant constants ---');
disp(['tau_c  = ', num2str(tau_c,  '%6.0f'), ' s']);
disp(['tau_r  = ', num2str(tau_r,  '%6.1f'), ' s']);
disp(['K_cool = ', num2str(K_cool, '%6.1f'), ' K per unit']);
disp(['K_load = ', num2str(K_load, '%6.2f'), ' K per %']);
disp(['y0     = ', num2str(y0,     '%6.2f')]);

%% PID gains from the SIMC tuning rule (Skogestad, 2003)
% For k/((tau1*s+1)(tau2*s+1)) with no time delay:
% Kc = tau1/(k*tau_cl),  tauI = min(tau1, 4*tau_cl),  tauD = tau2
tau_cl = 20;                            % chosen closed-loop time constant [s]
Kc   = tau_c/(K_cool*tau_cl);
tauI = min(tau_c, 4*tau_cl);
tauD = tau_r;
% Convert from SIMC (series) form to the parallel form Kp + Ki/s + Kd*s
Kp = Kc*(1 + tauD/tauI);
Ki = Kp/(tauI + tauD);
Kd = Kp*tauI*tauD/(tauI + tauD);
Tf = 1;                                 % derivative filter time constant [s]

disp(['Kp = ', num2str(Kp, '%6.3f'), '   Ki = ', num2str(Ki, '%7.4f'), '   Kd = ', num2str(Kd, '%5.2f')]);
disp(' ');

%% Define transfer functions of the plant
% Cooling path: command y -> radiator -> engine coolant
num = K_cool;
den = conv([tau_c 1],[tau_r 1]);   % (tau_c*s+1)(tau_r*s+1)
G = tf(num, den)

% Disturbance path: engine load -> engine coolant
numd = K_load;
dend = [tau_c 1];                  % (tau_c*s+1)
Gd = tf(numd, dend)

%% Define the three controllers
C_P   = tf(Kp, 1);                                    % Kp
C_PI  = tf([Kp Ki], [1 0]);                           % Kp + Ki/s
C_PID = tf(Kp, 1) + tf(Ki, [1 0]) + tf([Kd 0], [Tf 1]); % Kp + Ki/s + Kd*s/(Tf*s+1)

%% Closed loop transfer functions
% Error is e = T_c - T_ref, so a hotter engine gives more cooling.
% Cooling lowers T_c, so the two minus signs cancel and the
% standard negative feedback form can be used with G positive.
% Set-point -> T_c :  C*G / (1 + C*G)
T_P   = feedback(C_P*G, 1);
T_PI  = feedback(C_PI*G, 1);
T_PID = feedback(C_PID*G, 1);
% Engine load -> T_c :  Gd / (1 + C*G)
D_P   = Gd * feedback(1, C_P*G);
D_PI  = Gd * feedback(1, C_PI*G);
D_PID = Gd * feedback(1, C_PID*G);

%% Poles of the closed loop (stability check)
poles_P   = pole(T_P)
poles_PI  = pole(T_PI)
poles_PID = pole(T_PID)

%% Experiment 1 : set-point step 90 -> 95 °C
t = 0:1:600;               % time vector [s]
dT = 5;                    % size of set-point step [K]
[y_P,   t] = step(dT*T_P,   t);
[y_PI,  t] = step(dT*T_PI,  t);
[y_PID, t] = step(dT*T_PID, t);

figure;
plot(t, T_ref + y_P, t, T_ref + y_PI, t, T_ref + y_PID, 'LineWidth', 1.5);
hold on;
plot(t, (T_ref + dT)*ones(size(t)), 'k--');
grid on;
title('Set-point step 90 \rightarrow 95 °C');
xlabel('Time (s)');
ylabel('Coolant temperature T_c (°C)');
legend('P only', 'PI', 'PID', 'Set-point', 'Location', 'southeast');

% Step response characteristics
info_P   = stepinfo(T_P);
info_PI  = stepinfo(T_PI);
info_PID = stepinfo(T_PID);

disp('--- Experiment 1 : set-point step 90 -> 95 °C ---');
disp('Ctrl   Rise(s)  Peak(s)  Overshoot(%)  Settle(s)  Offset(K)');
disp(['P      ', num2str(info_P.RiseTime, '%6.1f'), '     ', num2str(info_P.PeakTime, '%6.1f'), '     ', ...
      num2str(info_P.Overshoot, '%6.1f'), '         ', num2str(info_P.SettlingTime, '%6.1f'), '      ', num2str(dT - y_P(end), '%5.2f')]);
disp(['PI     ', num2str(info_PI.RiseTime, '%6.1f'), '     ', num2str(info_PI.PeakTime, '%6.1f'), '     ', ...
      num2str(info_PI.Overshoot, '%6.1f'), '         ', num2str(info_PI.SettlingTime, '%6.1f'), '      ', num2str(dT - y_PI(end), '%5.2f')]);
disp(['PID    ', num2str(info_PID.RiseTime, '%6.1f'), '     ', num2str(info_PID.PeakTime, '%6.1f'), '     ', ...
      num2str(info_PID.Overshoot, '%6.1f'), '         ', num2str(info_PID.SettlingTime, '%6.1f'), '      ', num2str(dT - y_PID(end), '%5.2f')]);

%% Experiment 2 : engine load step 30 -> 80 % (truck hits a hill)
dL = 50;                   % size of load step [%]
[d_P,   t] = step(dL*D_P,   t);
[d_PI,  t] = step(dL*D_PI,  t);
[d_PID, t] = step(dL*D_PID, t);

figure;
plot(t, T_ref + d_P, t, T_ref + d_PI, t, T_ref + d_PID, 'LineWidth', 1.5);
hold on;
plot(t, T_ref*ones(size(t)), 'k--');
grid on;
title('Engine load step 30 \rightarrow 80 %');
xlabel('Time (s)');
ylabel('Coolant temperature T_c (°C)');
legend('P only', 'PI', 'PID', 'Set-point', 'Location', 'east');

% Peak rise and steady-state offset
peak_P   = max(d_P);    offset_P   = d_P(end);
peak_PI  = max(d_PI);   offset_PI  = d_PI(end);
peak_PID = max(d_PID);  offset_PID = d_PID(end);

% Settling time (±0.5 K band): go backwards from the end until
% the response leaves the band, same method as the time-domain notes
n = length(t);
s = n; while abs(d_P(s)   - d_P(end))   < 0.5, s = s - 1; end; settle_P   = t(s);
s = n; while abs(d_PI(s)  - d_PI(end))  < 0.5, s = s - 1; end; settle_PI  = t(s);
s = n; while abs(d_PID(s) - d_PID(end)) < 0.5, s = s - 1; end; settle_PID = t(s);

disp(' ');
disp('--- Experiment 2 : engine load step 30 -> 80 % ---');
disp('Ctrl   Peak(K)  Settle(s)  Offset(K)');
disp(['P      ', num2str(peak_P,   '%5.2f'), '     ', num2str(settle_P,   '%5.0f'), '      ', num2str(offset_P,   '%5.2f')]);
disp(['PI     ', num2str(peak_PI,  '%5.2f'), '     ', num2str(settle_PI,  '%5.0f'), '      ', num2str(offset_PI,  '%5.2f')]);
disp(['PID    ', num2str(peak_PID, '%5.2f'), '     ', num2str(settle_PID, '%5.0f'), '      ', num2str(offset_PID, '%5.2f')]);

% Cooling command needed by the PID (split range 0-3)
% e = T_c - T_ref = d_PID, so y = y0 + C_PID * e
U_PID = C_PID * D_PID;
[u_PID, t] = step(dL*U_PID, t);
disp(['PID cooling command: ', num2str(y0, '%4.2f'), ' -> ', num2str(y0 + u_PID(end), '%4.2f'), ...
      '  (max ', num2str(y0 + max(u_PID), '%4.2f'), ', limit 3)']);
disp('Split range: 0-1 valve, 1-2 pump, 2-3 fan');

% Inference:
% P only leaves a steady-state offset (~4 K for the load step), because
% it needs an error to produce any cooling. Adding I removes the offset
% (T_c returns to exactly 90 °C). Adding D lowers the peak and the
% overshoot. All closed-loop poles lie in the left half of the s-plane,
% so all three loops are stable.
