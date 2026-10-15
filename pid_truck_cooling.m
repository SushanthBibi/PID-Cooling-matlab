%% PID Thermal Control for Heavy-Duty Trucks
% Keeps the engine coolant temperature T_c at 90 °C.
% Transfer functions, step response, P vs PI vs PID
% (The Simulink model is built by hand from the same numbers)
% Sushanth Bibi (24BME0501), Ganishk Mohanraj (24BME0500)

clc;
clear;
close all;

%% Plant parameters (assumed example values)
% Each part is a first-order thermal lag K/(tau*s+1), the standard
% model of a thermal system (Ogata, Modern Control Engineering).
T_ref  = 90;     % set-point [°C]
tau_c  = 60;     % engine + coolant time constant [s]
tau_r  = 10;     % radiator time constant [s]
K_cool = 20;     % fall in T_c for one unit of cooling command [K]
K_load = 0.6;    % rise in T_c for 1 % more engine load [K/%]
y0     = 0.8;    % cooling command at 90 °C and 30 % load (scale 0-3)

%% PID gains from the SIMC tuning rule (Skogestad, 2003)
% For K/((tau1*s+1)(tau2*s+1)) with no time delay:
% Kc = tau1/(K*tau_cl),  tauI = min(tau1, 4*tau_cl),  tauD = tau2
tau_cl = tau_r;                         % closed-loop time constant [s]
Kc   = tau_c/(K_cool*tau_cl);           % = 0.3
tauI = min(tau_c, 4*tau_cl);            % = 40 s
tauD = tau_r;                           % = 10 s
% Convert from SIMC (series) form to parallel form Kp + Ki/s + Kd*s
Kp = Kc*(1 + tauD/tauI);                % = 0.375
Ki = Kp/(tauI + tauD);                  % = 0.0075
Kd = Kp*tauI*tauD/(tauI + tauD);        % = 3.0
Tf = 0.01;                              % derivative filter [s] (Simulink default N = 100)

%% Define transfer functions of the plant
% Cooling path: command y -> radiator -> engine coolant
num = K_cool;
den = conv([tau_c 1],[tau_r 1]);   % (60s+1)(10s+1)
G = tf(num, den)

% Disturbance path: engine load -> engine coolant
numd = K_load;
dend = [tau_c 1];                  % (60s+1)
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
% P only leaves a steady-state offset (~3.5 K for the load step), because
% it needs an error to produce any cooling. Adding I removes the offset
% (T_c returns to exactly 90 °C). Adding D lowers the peak and the
% overshoot. All closed-loop poles lie in the left half of the s-plane,
% so all three loops are stable.
