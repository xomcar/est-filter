clear vars
clc
close all

%% checkTFStability

fprintf("\n\t======\n");
z = tf('z');
G1 = (z^2-0.2)/(z^2-0.1*z);
G2 = (z + 2) / (z^2 + 0.4*z - 0.45);
G3 = (z^2 - 0.7*z + 1) / (z^2 + 0.4*z - 0.45);
G4 = z * (z + 1) / (z^2 + 2.5*z + 1);
Gs = [G1, G2, G3, G4];
bools = zeros(4,2); 

fprintf("Evaluating stability of functions:\n");

for i = 1:length(Gs)
    [a, b] = checkTFStability(Gs(i));
    bools(i,:) = [a, b];
    if([bools(i,:)] == [1, 1])
        fprintf("G%u is stable and has marginally stable inverse\n", i);
    elseif ([bools(i,:)] == [1, 0])
        fprintf("G%u is stable and has NOT marginally stable inverse\n", i);
    else
        fprintf("G%u is not stable\n", i);
    end
end

%% plotSpectrum

fprintf("\n\t======\n");
for i = 1:length(Gs)
   fprintf("Trying to plot spectrum of G%d\n", i);
   plotSpectrum(Gs(i)); 
end

%% residues integral

fprintf("\n\t======\n")
n1 = [1 0.8];
d1 = [1 0.3];
G11 = tf(n1, d1, -1);
n2 = [1 -0.8];
d2 = [1 -0.6 0.09];
G12 = tf(n2, d2, -1);
integ1 = resIntegral(G11);
integ2 = resIntegral(G12);
fprintf("Integral is %.2f, sholud be = 1\n", integ1);
fprintf("Integral is %.2f, should be = 0\n", integ2);
%% causal part

fprintf("\n\t======\n")
fprintf("Evaluating causal part and strictly anticausal part:\n");
n10 = [2 -1 -5];
d10 = [2 -7 3];
G10 = tf(n10, d10, -1);
[caus, sacaus] = causalPart(G10)
sum = minreal(caus + sacaus)
original = minreal(G10)