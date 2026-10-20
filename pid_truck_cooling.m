%% PID Temperature Control of a Truck Engine
% How the coolant temperature changes when the engine load goes up,
% first with no controller, then with a PID controller.
% Sushanth Bibi (24BME0501), Ganishk Mohanraj (24BME0500)

clc;
clear;
close all;

%% Plant (values from the references on the last slide)
% Fan path:  fan air flow [kg/s] -> coolant temperature [K]
num = 38.1;
den = conv([851 1],[11.3 1]);   % (851s+1)(11.3s+1)
G = tf(num, den)

% Load path: engine load [%] -> coolant temperature [K]
numd = 2.0;
dend = [851 1];                 % (851s+1)
Gd = tf(numd, dend)

%% PID controller (gains from the SIMC tuning rule)
Kp = 2.48;
Ki = 0.044;
Kd = 22.4;
C = tf(Kp, 1) + tf(Ki, [1 0]) + tf([Kd 0], [0.01 1]);   % Kp + Ki/s + Kd*s

%% Load goes from 30 % to 70 % at t = 0
t = 0:0.5:1800;                 % 30 minutes

% Without PID: load -> temperature = Gd
[T_noPID, t] = step(40*Gd, t);

% With PID: load -> temperature = Gd / (1 + C*G)
D = Gd * feedback(1, C*G);
[T_PID, t] = step(40*D, t);

%% Figure 1: no PID vs PID over 30 minutes
figure;
plot(t, 90 + T_noPID, 'r', t, 90 + T_PID, 'b', 'LineWidth', 2);
hold on;
plot(t, 90*ones(size(t)), 'k--');
plot(t, 105*ones(size(t)), 'm:', 'LineWidth', 1.5);
grid on;
title('Engine load 30 \rightarrow 70 %: no PID vs PID');
xlabel('Time (s)');
ylabel('Coolant temperature (°C)');
legend('No PID', 'With PID', 'Set-point 90 °C', 'Limit 105 °C', 'Location', 'northwest');

%% Figure 2: zoom on the moment the load changes
% Add 20 s of steady running at 90 °C before the load step (t = 0)
tz      = [-20; t(:)];
Tz_no   = 90 + [0; T_noPID(:)];
Tz_pid  = 90 + [0; T_PID(:)];

figure;
plot(tz, Tz_no, 'r', tz, Tz_pid, 'b', 'LineWidth', 2.5);
hold on;
plot(tz, 90*ones(size(tz)), 'k--');
plot([0 0], [89.8 93], 'g', 'LineWidth', 1.5);          % load step
grid on;
xlim([-20 150]);
ylim([89.8 93]);
title('Zoom: the moment the load goes up');
xlabel('Time (s)');
ylabel('Coolant temperature (°C)');
legend('No PID', 'With PID', 'Set-point 90 °C', 'Location', 'east');

%% Results
disp(['No PID:   temperature after 30 min = ', num2str(90 + T_noPID(end), '%5.1f'), ' °C']);
disp(['With PID: highest temperature      = ', num2str(90 + max(T_PID), '%5.2f'), ' °C']);
disp(['With PID: temperature after 30 min = ', num2str(90 + T_PID(end), '%5.2f'), ' °C']);
poles = pole(feedback(C*G, 1))   % all negative real parts -> stable

% Inference:
% Without a PID, the extra load keeps heating the engine: it passes the
% 105 °C limit after about 3 minutes and keeps rising. With the PID, the temperature rises by less
% than 1 °C and is brought back to exactly 90 °C in about 2 minutes.
