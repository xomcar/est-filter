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
    SigmaInf = A - GInf * C;
    
    %computing prediction 1 step
    X2Inf = A * X0 + GInf * (Y1 - C * X0);
    
    %computing error covariance 1 step
    P2Inf = SigmaInf * PInf * SigmaInf' + KInf * R * KInf' + QTilda; 
end

function [observable, detectable] = checkObsDetec(A, C, type)
    %%computing reachability
    ObsMat = obsv(A, C);
    if rank(ObsMat) == rank(A)
        observable = true;
        detectable = true;
        return
    else
        observable = false;
    end
    
    %%computing stabilizability
    %getting eigenvalues of A
    e = eig(A);
    
    %for each eigenvalue of A
    for i=1:length(e)
        %if continous system
        if type == "continuous"
            %check if eigenvalues is unstable
            if real(eig(i)) >= 0
                %evaluated PBH at eigen
                if rank(A) ~= rank([A - eig(i) * eye(size(A)); C])
                    detectable = false;
                    return;
                end
            end
        %if discrete system
        elseif type == "discrete"
            %check if eigenvalues is unstable
            if abs(eig(i)) >= 1
                %evaluated PBH at eigen
                if rank(A) ~= rank([A - eig(i) * eye(size(A)); C])
                    detectable = false;
                    return;
                end
            end
        end
    end
    
    detectable = true;
end

function [reachable, stabilizable] = checkReachStab(A, B, type)
    %%computing reachability
    ReachMat = ctrb(A, B);
    if rank(ReachMat) == rank(A)
        reachable = true;
        stabilizable = true;
        return
    else
        reachable = false;
    end
    
    %%computing stabilizability
    %getting eigenvalues of A
    e = eig(A);
    
    %for each eigenvalue of A
    for i=1:length(e)
        %if continous system
        if type == "continuous"
            %check if eigenvalues is unstable
            if real(eig(i)) >= 0
                %evaluated PBH at eigen
                if rank(A) ~= rank([A - eig(i) * eye(size(A)) B])
                    stabilizable = false;
                    return;
                end
            end
        %if discrete system
        elseif type == "discrete"
            %check if eigenvalues is unstable
            if abs(eig(i)) >= 1
                %evaluated PBH at eigen
                if rank(A) ~= rank([A - eig(i) * eye(size(A)) B])
                    stabilizable = false;
                    return;
                end
            end
        end
    end
    
    stabilizable = true;
end

