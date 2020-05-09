clearvars
clc

%% Parameters
MEASURES = 100;
NOISE_VARIANCE = 100;
NOISE_AVG = 0;
STATE_PARAMS = 2;
ZERO1 = - 0.9;
ZERO2 = -1.15;
POLE1 = 0.8 + 0.1j;
POLE2 = 0.8 - 0.1j;

%% Setting system and realization
W = tf(poly([ZERO1 ZERO2]), poly([POLE1 POLE2]), -1);
n = rand(1,MEASURES);
t = 0:MEASURES - 1;
x0 = NOISE_AVG + sqrt(NOISE_VARIANCE)*rand(STATE_PARAMS,1);
Sys = ss(W);
y = lsim(Sys, n, t, x0)';
[~, ~, C, ~] = ssdata(Sys);

%% Applying Filters
XKalman = zeros(STATE_PARAMS,MEASURES);
XKalmanSS = zeros(STATE_PARAMS,MEASURES);
YKalmanSS = zeros(STATE_PARAMS,MEASURES);
[XKalman(:, 1), P] = predKalman(Sys, y(1), x0, NOISE_VARIANCE * eye(STATE_PARAMS));
[XKalmanSS(:, 1), ~] = predKalmanSS(Sys, y(1), x0);
for i = 2:MEASURES
    [XKalman(:, i), P] = predKalman(Sys, y(i), XKalman(:, i - 1), P);
    [XKalmanSS(:, i), PSS] = predKalmanSS(Sys, y(i), XKalmanSS(:, i - 1));
end
[YWiener, varErrW] = WienerPredictor(W, y, 1);

%% Plotting
close all
figure
hold on
plot(t, y(1,:));
plot(t, C*XKalman);
plot(t, C*XKalmanSS);
plot(t, YWiener);
legend("Misure","Kalman","KalmanSS","Wiener")