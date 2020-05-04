function [estX, var_err] = distribMMSE(Y1, Y2, P, R1, R2, S1, S2)
%computing inverse befohand in order to save time
Pi = inv(P);
R1i = inv(R1);
R2i = inv(R2);

%computing these now in order to avoid computing them twice later
T1 = S1' * R1i * S1;
T2 = S2' * R2i * S2;

%same as X1 = centralMMSE(Y1, P, R1, S1);
var_err1 = inv(Pi + S1' * R1i * S1);
X1 = var_err1 * S1' * R1i * Y1;

%same as X2 = centralMMSE(Y2, P, R2, S2);
var_err2 = inv(Pi + S2' * R2i * S2);
X2 = var_err2 * S2' * R2i * Y2;

%implementing formulas
var_err = inv(Pi + T1 + T2);
estX = var_err * (T1 * X1 + T2 * X2);
end