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
