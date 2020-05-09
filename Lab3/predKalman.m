function [X2, P2] = predKalman(sys, Y1, X0, P0)
    [A, B, C, D] = ssdata(sys);
    
    %supposing unit variance white noise we get a noise covariance:
    Q = B*B';
    R = D*D';
    S = B*D';
    
    %after uncorrelation of noises
    QTilda = Q - S / R * S';
    F = A - S / R * C;
    
    %computing various time variant matrices
    L1 = P0 * C' / (C * P0 * C' + R);
    K1 = F * L1;
    G1 = K1 + S / R;
    Sigma1 = A - G1 * C;
    
    %computing prediction 1 step
    X2 = A * X0 + G1 * (Y1 - C * X0);
    
    %computing error covariance 1 step
    P2 = Sigma1 * P0 * Sigma1' + K1 * R * K1' + QTilda;
end