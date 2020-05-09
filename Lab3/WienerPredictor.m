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