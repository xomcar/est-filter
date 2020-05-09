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