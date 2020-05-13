clearvars
clc

%% Parameters
MEASURES = 100;
MEASURES_PARAMS = 1;
NOISE_VARIANCE = 100;
NOISE_AVG = 0;
STATE_PARAMS = 2;
ZERO1 = - 0.9;
ZERO2 = -1.15;
POLE1 = 0.8 + 0.1j;
POLE2 = 0.8 - 0.1j;
ERR_VAR0 = 1;

%% Setting system and realization
W = tf(poly([ZERO1 ZERO2]), poly([POLE1 POLE2]), -1);
n = randn(MEASURES_PARAMS,MEASURES);
t = 0:MEASURES - 1;
x0 = NOISE_AVG + sqrt(NOISE_VARIANCE)*randn(STATE_PARAMS,1);
Sys = ss(W);
y = lsim(Sys, n, t, x0)';
[~, ~, C, D] = ssdata(Sys);

%% Applying Filters
XKalman = zeros(STATE_PARAMS,MEASURES);
XKalmanSS = zeros(STATE_PARAMS,MEASURES);
YKalmanSS = zeros(MEASURES_PARAMS,MEASURES);
[XKalman(:, 1), P] = predKalman(Sys, y(1), x0, ERR_VAR0 * eye(STATE_PARAMS));
[XKalmanSS(:, 1), ~] = predKalmanSS(Sys, y(1), x0);
for i = 2:MEASURES
    [XKalman(:, i), P] = predKalman(Sys, y(i), XKalman(:, i - 1), P);
    [XKalmanSS(:, i), PSS] = predKalmanSS(Sys, y(i), XKalmanSS(:, i - 1));
end
[YWiener, varErrW] = WienerPredictor(W*W', y, 1);
YWiener = real(YWiener');

%% Plotting
close all
figure()
subplot(3,1,1)
hold on
plot(t, y(1,:), 'r');
plot(t, C*XKalman, 'g');
legend("Realization", "Kalman 1 step")
subplot(3,1,2)
hold on
plot(t, y(1,:), 'r');
plot(t, C*XKalmanSS, 'b');
legend("Realization", "KalmanSS 1 step")
subplot(3,1,3)
hold on
plot(t, y(1,:), 'r');
plot(t, YWiener, 'k');
legend("Realization", "Wiener 1 step")

%% Comparing Wiener as a particular KalmanSS case
R = D*D';
LambdaInf = C * PSS * C'+ R;
LamdaInfSQRT = sqrt(LambdaInf);

%% Computing and printing transfer functions
WN = tf(Sys) * LamdaInfSQRT
L = spectralFactor(W*W')

%% Plotting
figure()
hold on
plot(impulse(L));
plot(impulse(WN));
legend("MP spectral factor of W*W'", "Transfer function of innovation on KSS")


%% Local function

function [pred,varErr] = WienerPredictor(tranFun, realizY, kStep)

%computing spectral factor
L = spectralFactor(tranFun);
%L = spectral(tranFun);
%extracting num and denum from L
[C, A] = tfdata(L, 'v');
%polynomial division
%not needed to convert poly into roots because
%we are adding k 0*z^i items to the right
%and parameters are in decreasing order
[Q, R] = deconv([C, zeros(1, kStep)], [A, 0]);

%remainder gives Ck(z), Lk = Ck / A
Lk = minreal(tf(R, A, -1));

%quotient gives anticausal part
%from anticausal coefficients compute variance
varErr = 0;
for k = 1:length(Q)
    varErr = varErr + (Q(k))^2;
end
t = 0:length(realizY)-1;

%filter system with after 1/L 
H = Lk / L;
pred = lsim(H, realizY, t);
end

function L = spectralFactor(tranFun)
%extract data from zero pole gain
[zeros, poles, gain] = zpkdata(tranFun, 'v');

%remove '0' zeros
zeros = nonzeros(zeros);

%initial values
lambda2 = gain;
L = 1;

%numerator denumerator degrees
m = 0; n = 0;
z = tf('z');

%extract marg stab zeroes, compute numerator of lambda
for i = 1:length(zeros)
    %if zero is marginally stable
    if abs(zeros(i)) <= 1
        %add zero to transfer function
        L = L * (z - zeros(i));
        %compute zero in numerator count
        m = m + 1;
        %update lambda2
        lambda2 = lambda2 * (-zeros(i));
    end
end

%extract stab poles, compute denumerator of lambda
for i = 1:length(poles)
    %if pole is stable
    if abs(poles(i)) < 1 
        %add poles to transfer function
        L = L /(z - poles(i));
        %compute pole in denumerator count
        n = n + 1;
        %update lambda2
        lambda2 = lambda2 / (-poles(i));
    end
end

%compute lambda
lambda = sqrt(lambda2);

%compute spectral factorization
L = z^(n-m) * lambda * L;
end

