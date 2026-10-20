%% PID Temperature Control of a Truck Engine
% How the coolant temperature changes when the engine load goes up,
% first with no controller, then with a PID controller.
% Sushanth Bibi (24BME0501), Ganishk Mohanraj (24BME0500)
%
% Sources
% [1] Volvo Trucks, Fact sheet D13C460 EU5SCR-M, 2010
% [2] G. S. Rosa et al., "Dimensioning procedure of coolant radiators for
%     trucks and buses", COBEM 2009
% [3] M. Mohamed et al., "Performance of electronically controlled
%     automotive engine cooling system using PID and LQR control
%     techniques", IOSR J. Engineering 8(6), 2018, pp. 16-25
% [4] A. Haury, J. Volkering, "Modelisation of the engine coolant
%     warming-up behavior", MSc thesis 2011:13, Chalmers
% [5] Engineering ToolBox, specific heat and density tables
% [6] S. Skogestad, "Simple analytic rules for model reduction and PID
%     controller tuning", J. Process Control 13 (2003) 291-309

clc;
clear;
close all;

%% Data from the sources
P     = 338;        % engine power [kW]                       [1]
m_eng = 1121;       % engine dry weight [kg]                  [1]
V_c   = 38;         % coolant volume [L]                      [1]
V_o   = 33;         % oil volume [L]                          [1]
cp_fe = 0.46;       % cast iron [kJ/kg K]                     [5]
rho_w = 0.96531;    % water at 90 °C [kg/L]                   [5]
cp_w  = 4.2053;     % water at 90 °C [kJ/kg K]                [5]
rho_o = 0.900;      % lubricating oil [kg/L]                  [5]
cp_o  = 1.67;       % mineral oil [kJ/kg K]                   [5]
cp_a  = 1.006;      % air [kJ/kg K]                           [5]
beta  = 0.5;        % heat to coolant / engine power (0.45-0.55) [2]
m_w   = 3.25;       % coolant flow [kg/s]                     [2]
T_w   = 100;        % radiator rating: water in [°C]          [2]
T_a   = 30;         % radiator rating: air in [°C]            [2]
% Fan air flow = where system curve meets fan curve, Fig. 6 of [2]:
%   system: 25.105x^2 + 76.084x - 3.7131,   fan: -1.9147x^2 - 457.51x + 2945
mf_max = max(roots([25.105+1.9147, 76.084+457.51, -3.7131-2945]));   % = 4.50 kg/s

%% Transfer function numbers from the energy balance [3]
% Block + coolant + oil as one mass [4]:
%   C*dT/dt = Q_engine - eps*m_fan*cp_a*(T - T_a)
C     = m_eng*cp_fe + V_c*rho_w*cp_w + V_o*rho_o*cp_o;   % heat capacity [kJ/K]
eps   = beta*P/(mf_max*cp_a*(T_w - T_a));               % radiator effectiveness
H     = beta*P*0.30/(90 - T_a);                         % heat removed per K at 90 °C, 30 % load [kW/K]
tau_c = C/H;                                            % engine time constant [s]
tau_r = V_c*rho_w/m_w;                                  % coolant loop time [s]
K_fan = eps*cp_a*(90 - T_a)/H;                          % K per kg/s of fan air
K_load = beta*P/100/H;                                  % K per % load

%% Plant
% Fan path:  fan air flow [kg/s] -> coolant temperature [K]
num = K_fan;                        % = 38.1
den = conv([tau_c 1],[tau_r 1]);    % (851s+1)(11.3s+1)
G = tf(num, den)

% Load path: engine load [%] -> coolant temperature [K]
numd = K_load;                      % = 2.0
dend = [tau_c 1];                   % (851s+1)
Gd = tf(numd, dend)

%% PID controller (SIMC tuning rule [6])
Kc   = tau_c/(K_fan*tau_r);
tauI = min(tau_c, 4*tau_r);
tauD = tau_r;
Kp = Kc*(1 + tauD/tauI)             % = 2.48
Ki = Kp/(tauI + tauD)               % = 0.044
Kd = Kp*tauI*tauD/(tauI + tauD)     % = 22.4
C = tf(Kp, 1) + tf(Ki, [1 0]) + tf([Kd 0], [0.01 1]);   % Kp + Ki/s + Kd*s

%% Engine load: 30 % for 20 s, then rises linearly to 70 % over 60 s
t = (0:0.5:1800)';              % 30 minutes
t_ramp = 60;                    % time for the load to rise [s]
dL = 40*min(max(t - 20, 0)/t_ramp, 1);   % change in load from 30 % [%]
load_pct = 30 + dL;             % engine load [%]

% Without PID: load -> temperature = Gd
T_noPID = lsim(Gd, dL, t);

% With PID: load -> temperature = Gd / (1 + C*G)
D = Gd * feedback(1, C*G);
T_PID = lsim(D, dL, t);

%% Figure 1: 30 minutes
figure;
yyaxis left;
plot(t, 90 + T_noPID, 'r-', t, 90 + T_PID, 'b-', 'LineWidth', 2);
hold on;
plot(t, 90*ones(size(t)), 'k--', t, 105*ones(size(t)), 'm:');
ylabel('Coolant temperature (°C)');
ylim([80 180]);
yyaxis right;
plot(t, load_pct, 'g-', 'LineWidth', 1.5);
ylabel('Engine load (%)');
ylim([0 100]);
grid on;
xlim([0 1800]);
title('No PID vs PID');
xlabel('Time (s)');
legend('No PID', 'With PID', 'Set-point 90 °C', 'Limit 105 °C', 'Engine load', 'Location', 'east');

%% Figure 2: zoom on the load increase
figure;
yyaxis left;
plot(t, 90 + T_noPID, 'r-', t, 90 + T_PID, 'b-', 'LineWidth', 2.5);
hold on;
plot(t, 90*ones(size(t)), 'k--');
ylabel('Coolant temperature (°C)');
ylim([89.8 93]);
yyaxis right;
plot(t, load_pct, 'g-', 'LineWidth', 1.5);
ylabel('Engine load (%)');
ylim([0 100]);
grid on;
xlim([0 200]);
title('Zoom: while the load goes up');
xlabel('Time (s)');
legend('No PID', 'With PID', 'Set-point 90 °C', 'Engine load', 'Location', 'east');

%% Results
disp(['No PID:   temperature after 30 min = ', num2str(90 + T_noPID(end), '%5.1f'), ' °C']);
disp(['With PID: highest temperature      = ', num2str(90 + max(T_PID), '%5.2f'), ' °C']);
disp(['With PID: temperature after 30 min = ', num2str(90 + T_PID(end), '%5.2f'), ' °C']);
poles = pole(feedback(C*G, 1))   % all negative real parts -> stable

% Inference:
% Without a PID, the extra load keeps heating the engine: it passes the
% 105 °C limit after about 4 minutes and keeps rising. With the PID, the
% temperature rises by only 0.6 °C while the load goes up, and is back at
% exactly 90 °C about 80 s after the load stops rising.
