function integ = resIntegral(tranFun)
%returns integral of tranFun computed on the unit circle via residue method

z = tf('z');
tranFun = tranFun / z;
%extract numerator and denumerator from transfer function obj
[num, den] = tfdata(tranFun, 'v');
%determine residues and decomposition
[res, poles] = residue(num,den);
%res'
%poles'
integ = 0;
i = 2;
while(i <= length(poles) + 1)
    %if pole is stable
    if abs(poles(i-1)) < 1
        %add to integral corresponding residue
        integ = integ + res(i-1);
        %fprintf("Added %f from pole %f\n", res(i-1), poles(i-1));
    end
    %skip multiplicity
    while i <= length(poles) && poles(i-1) == poles(i)
        i = i + 1;
        %fprintf("Skipped %f from pole %f\n", res(i-1), poles(i-1));
    end
    i = i + 1;
end
end