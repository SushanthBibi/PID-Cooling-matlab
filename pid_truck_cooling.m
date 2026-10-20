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

%% 1) No controller: load goes from 30 % to 70 %
t1 = 0:1:3600;                  % one hour
[T_open, t1] = step(40*Gd, t1);

figure;
plot(t1, 90 + T_open, 'LineWidth', 1.5);
grid on;
title('No controller: engine load 30 \rightarrow 70 %');
xlabel('Time (s)');
ylabel('Coolant temperature (°C)');

%% 2) With PID: same load step
% Closed loop: load -> temperature = Gd / (1 + C*G)
D = Gd * feedback(1, C*G);
t2 = 0:0.1:300;                 % five minutes
[T_pid, t2] = step(40*D, t2);

figure;
plot(t2, 90 + T_pid, 'LineWidth', 1.5);
hold on;
plot(t2, 90*ones(size(t2)), 'k--');
grid on;
title('With PID: engine load 30 \rightarrow 70 %');
xlabel('Time (s)');
ylabel('Coolant temperature (°C)');
legend('Coolant temperature', 'Set-point 90 °C');

%% Results
disp(['No controller: temperature after 1 hour = ', num2str(90 + T_open(end), '%5.1f'), ' °C']);
disp(['With PID: highest temperature = ', num2str(90 + max(T_pid), '%5.2f'), ' °C']);
disp(['With PID: final temperature   = ', num2str(90 + T_pid(end), '%5.2f'), ' °C']);
poles = pole(feedback(C*G, 1))   % all negative real parts -> stable

% Inference:
% Without a controller, the extra load keeps heating the engine until
% it is far above 90 °C. With the PID, the temperature rises by less
% than 1 °C and is brought back to exactly 90 °C in about 2 minutes.
