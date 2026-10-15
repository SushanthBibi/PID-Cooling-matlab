%% PID Thermal Control for Heavy-Duty Trucks
% Keeps the engine coolant temperature T_c at 90 °C.
% Transfer functions, step response, P vs PI vs PID
% (The Simulink model is built by hand from the same numbers)
% Sushanth Bibi (24BME0501), Ganishk Mohanraj (24BME0500)
%
% References (numbers in square brackets below)
% [1] Volvo Trucks, Fact sheet D13C460 EU5SCR-M, 2010 (engine data)
% [2] G. S. Rosa et al., "Dimensioning procedure of coolant radiators for
%     trucks and buses", COBEM 2009 (heat rejection, radiator, fan)
% [3] M. Mohamed et al., "Performance of electronically controlled
%     automotive engine cooling system using PID and LQR control
%     techniques", IOSR J. Engineering 8(6), 2018, pp. 16-25 (model)
% [4] A. Haury, J. Volkering, "Modelisation of the engine coolant
%     warming-up behavior", MSc thesis 2011:13, Chalmers (lumped mass)
% [5] Engineering ToolBox, property tables (specific heat, density)
% [6] S. Skogestad, "Simple analytic rules for model reduction and PID
%     controller tuning", J. Process Control 13 (2003) 291-309 (SIMC)
% [7] N. S. Nise, Control Systems Engineering (feedback, final value thm)
% [8] MATLAB/Simulink, PID Controller block (default filter N = 100)

clc;
clear;
close all;

%% Data from the references
% Volvo D13C460 truck engine [1]
P_rated = 338;      % rated power [kW]
m_eng   = 1121;     % dry weight [kg]
V_cool  = 38;       % cooling system volume [L]
V_oil   = 33;       % oil volume [L]
% Material properties [5]
cp_iron = 0.46;     % cast iron [kJ/kg K]
rho_w   = 0.96531;  % water at 90 °C [kg/L]
cp_w    = 4.2053;   % water at 90 °C [kJ/kg K]
rho_oil = 0.900;    % lubricating oil [kg/L]
cp_oil  = 1.67;     % mineral oil [kJ/kg K]
cp_air  = 1.006;    % air [kJ/kg K]
% Truck/bus cooling system [2]
beta    = 0.50;     % heat to coolant / engine power (Eq. 6: 0.45-0.55, middle)
m_cool  = 3.25;     % coolant mass flow through the radiator [kg/s]
T_wr    = 100;      % radiator rating: water inlet [°C] (373.15 K)
T_air   = 30;       % radiator rating: air inlet [°C] (303.15 K)
% Fan operating point = system curve meets fan curve (Fig. 6 of [2])
r      = roots([25.105 + 1.9147, 76.084 + 457.51, -3.7131 - 2945]);
mf_max = max(r);    % largest fan air flow [kg/s] (= 4.50)

%% Test conditions
T_ref = 90;         % set-point [°C] (project requirement, slide 2)
L0    = 30;         % engine load before the step [%]

%% Plant constants from the energy balance [3]
% Engine, coolant and oil lumped as one mass [4]:
%   C*dT/dt = Q_in - eps*mf*cp_air*(T - T_air)       (Eqs. 1-3 of [3])
C   = m_eng*cp_iron + V_cool*rho_w*cp_w + V_oil*rho_oil*cp_oil;   % [kJ/K]
Q_L = beta*P_rated/100;                  % heat to coolant per % load [kW/%]
eps = beta*P_rated/(mf_max*cp_air*(T_wr - T_air));  % radiator effectiveness
mf0 = Q_L*L0/(eps*cp_air*(T_ref - T_air));          % fan flow at 90 °C, 30 % [kg/s]
H0  = eps*cp_air*mf0;                    % heat removed per K [kW/K]

tau_c  = C/H0;                           % engine time constant [s]
K_cool = eps*cp_air*(T_ref - T_air)/H0;  % fall in T_c per kg/s of fan air [K]
K_load = Q_L/H0;                         % rise in T_c per % load [K/%]
% Coolant needs time to go round the loop: delay = coolant mass / flow.
% It is written as a first-order lag 1/(tau_r*s+1)
tau_r  = V_cool*rho_w/m_cool;            % loop time [s]
y0     = mf0;                            % fan air flow at the operating point [kg/s]
y_max  = mf_max;                         % fan limit [kg/s]

disp('--- Plant constants (calculated) ---');
disp(['C      = ', num2str(C,      '%6.1f'), ' kJ/K']);
disp(['tau_c  = ', num2str(tau_c,  '%6.0f'), ' s']);
disp(['tau_r  = ', num2str(tau_r,  '%6.2f'), ' s']);
disp(['K_cool = ', num2str(K_cool, '%6.2f'), ' K per kg/s']);
disp(['K_load = ', num2str(K_load, '%6.3f'), ' K per %']);
disp(['y0     = ', num2str(y0,     '%6.3f'), ' kg/s (limit ', num2str(y_max, '%4.2f'), ')']);

%% PID gains from the SIMC tuning rule [6]
% For K/((tau1*s+1)(tau2*s+1)):
% Kc = tau1/(K*tau_cl),  tauI = min(tau1, 4*tau_cl),  tauD = tau2
tau_cl = tau_r;                         % tau_cl = delay, "tight control" in [6]
Kc   = tau_c/(K_cool*tau_cl);
tauI = min(tau_c, 4*tau_cl);
tauD = tau_r;
% Convert from SIMC (series) form to parallel form Kp + Ki/s + Kd*s [6]
Kp = Kc*(1 + tauD/tauI);
Ki = Kp/(tauI + tauD);
Kd = Kp*tauI*tauD/(tauI + tauD);
Tf = 0.01;                              % derivative filter [s], Tf = 1/N, N = 100 [8]
disp(['Kp = ', num2str(Kp, '%6.3f'), '   Ki = ', num2str(Ki, '%7.4f'), '   Kd = ', num2str(Kd, '%6.2f')]);
disp(' ');

%% Define transfer functions of the plant
% Cooling path: fan air flow -> radiator -> engine coolant
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
% feedback(A,B) = A/(1 + A*B), negative feedback [3]
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

%% Experiment 2 : engine load step 30 -> 70 % (truck hits a hill)
dL = 40;                   % size of load step [%]
[d_P,   t] = step(dL*D_P,   t);
[d_PI,  t] = step(dL*D_PI,  t);
[d_PID, t] = step(dL*D_PID, t);

figure;
plot(t, T_ref + d_P, t, T_ref + d_PI, t, T_ref + d_PID, 'LineWidth', 1.5);
hold on;
plot(t, T_ref*ones(size(t)), 'k--');
grid on;
title('Engine load step 30 \rightarrow 70 %');
xlabel('Time (s)');
ylabel('Coolant temperature T_c (°C)');
legend('P only', 'PI', 'PID', 'Set-point', 'Location', 'east');

% Peak rise and steady-state offset
peak_P   = max(d_P);    offset_P   = d_P(end);
peak_PI  = max(d_PI);   offset_PI  = d_PI(end);
peak_PID = max(d_PID);  offset_PID = d_PID(end);

% Settling time (±0.1 K band): go backwards from the end until
% the response leaves the band, same method as the time-domain notes
n = length(t);
s = n; while abs(d_P(s)   - d_P(end))   < 0.1, s = s - 1; end; settle_P   = t(s);
s = n; while abs(d_PI(s)  - d_PI(end))  < 0.1, s = s - 1; end; settle_PI  = t(s);
s = n; while abs(d_PID(s) - d_PID(end)) < 0.1, s = s - 1; end; settle_PID = t(s);

disp(' ');
disp('--- Experiment 2 : engine load step 30 -> 70 % ---');
disp('Ctrl   Peak(K)  Settle(s)  Offset(K)');
disp(['P      ', num2str(peak_P,   '%5.2f'), '     ', num2str(settle_P,   '%5.0f'), '      ', num2str(offset_P,   '%5.2f')]);
disp(['PI     ', num2str(peak_PI,  '%5.2f'), '     ', num2str(settle_PI,  '%5.0f'), '      ', num2str(offset_PI,  '%5.2f')]);
disp(['PID    ', num2str(peak_PID, '%5.2f'), '     ', num2str(settle_PID, '%5.0f'), '      ', num2str(offset_PID, '%5.2f')]);

% Fan air flow needed by the PID
% e = T_c - T_ref = d_PID, so y = y0 + C_PID * e
U_PID = C_PID * D_PID;
[u_PID, t] = step(dL*U_PID, t);
disp(['PID fan air flow: ', num2str(y0, '%4.2f'), ' -> ', num2str(y0 + u_PID(end), '%4.2f'), ...
      ' kg/s  (max ', num2str(y0 + max(u_PID), '%4.2f'), ', limit ', num2str(y_max, '%4.2f'), ')']);

% Inference:
% P only leaves a steady-state offset, because it needs an error to
% produce any cooling. Final value theorem [7] (s = 0) for the load step:
% offset = dL*K_load/(1 + Kp*K_cool), about 0.8 K.
% Adding I makes C(0) infinite, so the offset becomes 0 (T_c returns to
% exactly 90 °C). Adding D lowers the load peak and the set-point
% overshoot. All closed-loop poles lie in the left half of the s-plane,
% so all three loops are stable.
