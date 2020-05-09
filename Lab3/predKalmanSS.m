function [X2Inf, P2Inf] = predKalmanSS(sys, Y1, X0)
    [A, B, C, D] = ssdata(sys);
    
    %supposing unit variance white noise we get a noise covariance:
    Q = B*B';
    R = D*D';
    S = B*D';
    
    %after uncorrelation of noises
    QTilda = Q - S / R * S';
    F = A - S / R * C;
    
     %check for validity of FK theorem
    [~, detec] = checkObsDetec(F, C, "discrete");
    [~, stab] = checkReachStab(F, sqrt(QTilda), "discrete");
    
    if ~(detec == 1 && stab == 1)
        X2Inf = [];
        P2Inf = [];
        return
    end
    
    %solve ARE
    PInf = dare(F', C', QTilda, R);
    KInf = F * PInf * C' / (C * PInf * C' + R);
    GInf = KInf + S / R;
    SigmaInf = F - KInf * C;
    
    %computing prediction 1 step
    X2Inf = A * X0 + GInf * (Y1 - C * X0);
    
    %computing error covariance 1 step
    P2Inf = SigmaInf * PInf * SigmaInf' + KInf * R * KInf' + QTilda;    
end