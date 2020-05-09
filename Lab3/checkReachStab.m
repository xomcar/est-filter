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

