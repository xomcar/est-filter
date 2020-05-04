 function [tfC, tfSAC] = causalPart(tranFun)
    %compute causal part and STRICTLY anticausal part of tranFun
    
    %% setting up an initial stable - unstable division
    %extract numerator and denumerator from transfer function obj
    [num, den] = tfdata(tranFun, 'v');
    
    %determine decomposition
    [res, poles, k] = residue(num,den);
    
    %preallocating matrix of poles-residues
    npols = length(poles);
    stableRes = zeros(npols, 1);
    stablePoles = zeros(npols, 1);
    unstableRes = zeros(npols, 1);
    unstablePoles = zeros(npols, 1);
    
    %dividing the poles by stability
    st = 1; un = 1;
    for i=1:length(poles)
        %if stable pole, add to 
        if abs(poles(i)) < 1
            stablePoles(st) = poles(i);
            stableRes(st) = res(i);
            st = st + 1;
        else
            unstablePoles(un) = poles(i);
            unstableRes(un) = res(i);
            un = un + 1;
        end
    end
    
    %vector cropping
    stableRes = stableRes(1:st-1);
    unstableRes = unstableRes(1:un-1);
    stablePoles = stablePoles(1:st-1);
    unstablePoles = unstablePoles(1:un-1);
    
    %% computing k0-
    %definining k0-, index variable and number of unstable poles
    k0m = 0; i = 1; nunst = length(unstablePoles);
    %for every unstable pole
    while i <= length(unstablePoles)
        %reset multiplicity
        u = 1;
        %evaluate each pole with its multiplicity
        while i <= nunst
            pole = unstablePoles(i);
            res = unstableRes(i);
            %compute at z = 0 and add to k0-
            k0m = k0m - (res / (0 - pole)^u);
            %if it has multiplicity, raise and continue
            if i ~= nunst && pole == unstablePoles(i + 1)
                i = i + 1;
                u = u + 1;
                continue
            else
                %if there are no more multiplicity, next new pole
                i = i + 1;
                break
            end
        end
    end
    
    %% computing k0+
    %if k is empty, k0 takes a 0 value
    k0 = 0;
    %else
    if ~isempty(k)
        %extract k0
        k0 = k(end);
    end
    %k0+ = k0 - k0-
    k0p = k0 - k0m;
    
    %% computing causal part
    [causalNum, causalDen] = residue(stableRes, stablePoles, k0p);
    tfC = tf(causalNum, causalDen, -1);
    
    %% computing strictly anticausal part
    tfSAC = minreal(tranFun - tfC);
end