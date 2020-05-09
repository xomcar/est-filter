function [internal, external] = checkStability(A, B, C, D, type)
       
    %%checking eigenvalues for internal stability
    
    %compute eigenvalues of a
    e = eig(A);
    for i = 1:length(e)
        if type == "continuous"
            %if system is continous
            if real(e(i)) >= 0
                %if internally unstable
                internal = false;
                external = false;
                return;
            end
        elseif type == "discrete"
            %if system is discrete
            if abs(e(i)) >= 1
                %if internally unstable
                internal = false;
                external = false;
                return;
            end
        else
            %impossible configuration used as error for input string
            internal = false;
            external = true;
            return;
        end
    end
    
    internal = true;
    
    %%checking poles for external stability
    
    %generating system depending on the type
    if type == "continuous"
        W = tf(ss(A, B, C, D));
    else
        W = tf(ss(A, B, C, D, -1));
    end
    
    %generating minimal realization and extracting poles
    W = minreal(W);
    [~, d] = tfdata(W);
    poles = roots(d);
    
    %evaluating poles for external stability
    for i = 1:length(poles)
        if type == "continuous"
            %if system is continous
            if real(poles(i)) >= 0
                %if externally unstable
                external = false;
                return;
            end
        else
            %if system is discrete
            if abs(poles(i)) >= 1
                %if externally unstable
                external = false;
                return;
            end
        end
    end
    
    external = true;
end

