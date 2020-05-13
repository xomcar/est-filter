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

