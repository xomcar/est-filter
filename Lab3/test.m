clearvars
clc

%% Setting system and realization
z = tf('z');
W = (z + 0.9)*(z + 1.15) / ((z - 0.8 + 0.1j)*(z - 0.8 - 0.1j));
n = rand(1,100);
t = 0:99;
x0 = sqrt(100)*rand(2,1);
Sys = ss(W);
[y, ~, x] = lsim(Sys, n, t, x0);
y = y'; x = x';
[A, B, C, D] = ssdata(Sys);

%% Applying Filters
XKalman = zeros(2,100);
XKalmanSS = zeros(2,100);
[XKalman(:, 1), P] = predKalman(Sys, y(1), x0, sqrt(100) * eye (2));
[XKalmanSS(:, 1), ~] = predKalmanSS(Sys, y(1), x0);
for i = 2:100
    [XKalman(:, i), P] = predKalman(Sys, y(i), XKalman(:, i - 1), P);
    [XKalmanSS(:, i), PSS] = predKalmanSS(Sys, y(i), XKalmanSS(:, i - 1));
end


%% Plotting
close all
figure
hold on
plot(t, y(1,:));
plot(t, C*XKalman);
plot(t, C*XKalmanSS);
legend("Misure","Kalman","KalmanSS")