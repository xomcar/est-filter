function [boolStable, boolMStable] = checkTFStability(tranFun)
    %compute stability and marginal stability of inverse of a tranFun
    
    %set output variables default as false
    boolStable = false;
    boolMStable = false;
    
    %extract num and denum from tranfun
    [N, D] = tfdata(tranFun, 'v');
    
    %find roots and zeroes
    zeroes = roots(N);
    %fprintf("Zeroes: %i\n",abs(zeroes))
    poles = roots(D);
    %fprintf("Poles: %i\n",abs(poles))
    
    %check stability of poles
    if (abs(poles) < 1)
        boolStable = true;
        %check marginal stability of zeroes
        if (abs(zeroes) <= 1)
            boolMStable = true;
        end
    end
end