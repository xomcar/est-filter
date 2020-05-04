function [estX, var_err] = centralMMSE(Y, P, R, S)
%computing inverse befohand in order to save time
Pi = inv(P);
Ri = inv(R);

%implementing formulas
var_err = inv(Pi + S' * Ri * S);
estX = var_err * S' * Ri * Y;
end