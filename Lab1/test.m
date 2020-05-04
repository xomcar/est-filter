%% part 1

clc
clear all

v = ones(1,10);
v = zeroTail(v);
fprintf("Zero'd last 3 item of 10-dim ones vec: \n\n\t[");
fprintf(' %g ', v);
fprintf("]\n\n");

fprintf('\t========================\n\n')

fprintf("Flipped diagonal on a 4x4 random matrix:\n\n");
mat = round(rand(4,4)*25);
disp(mat);
mat = flipDiag(mat);
disp(mat);

fprintf('\t========================\n\n')

fprintf("Created polinomial p(x) with roots {-1/3, 3/4 + j, 3/4 - j}:\n\n");
l = [-1/3, 3/4 + 1i, 3/4 - 1i];
p = poly(l);
fprintf("\t(%3.5f)x^3 + (%3.5f)x^2 + (%3.5f)x + (%3.5f)\n\n", p);
q = [1 -0.5 0];
fprintf("g(x) is: \n\n\t(%3.5f)x^2 + (%3.5f)x + (%3.5f)\n\n", q);
if (testSchur(p, q))
    fprintf("Polynomial p*g(x) is Schur stable.\n");
else
    fprintf("Polynomial p*g(x) is not Schur stable.\n");
end
fprintf("Plotting product in the [-10, 10] interval...\n\n");

%% part 2

clear all

meas = 500;
vars = 2;

fprintf('\t========================\n\n')
fprintf("Computing linear estimator for %d gaussian events \n", vars);
fprintf("\ton %d measures with gaussian noise...\n\n", meas);

%prefilters values
S1 = rand();
S1 = repmat(S1, meas, vars);
S2 = rand();
S2 = repmat(S2, meas, vars);
S = [S1; S2];

%priori variance as pos semidef simmetric
P = rand(vars);
P = P * P';

%noise variance matrices as pos semidef simmetric
R1 = rand(meas) * 0.01;
R1 = R1 * R1';
R2 = rand(meas) * 0.01;
R2 = R2 * R2';
R = blkdiag(R1, R2);

%noises are a 1000x2 random vector with zero mean and variance R
Noise_mean = zeros(vars * meas, 1);
Noise = mvnrnd(Noise_mean, R)';

%random event measured is a 2x1 random vector with zero mean and variance P
X_mean = zeros(vars,1);
X = mvnrnd(X_mean, P)';

%realizing Y as a 1000x1 random vector adding noises to filtered events
Y = S * X + Noise;
Y1 = Y(1:meas,:);
Y2 = Y(meas + 1:end,:);

%estimate computation
%central
tic
[x_est_c, covar_err_c] = centralMMSE(Y, P, R, S);
tc = toc;
fprintf('Time elapsed for cumulative MMSE: %.2f ms.\n', tc * 1000);

%distributed
tic
[x_est_d, covar_err_d] = distribMMSE(Y1, Y2, P, R1, R2, S1, S2);
td = toc;
fprintf('Time elapsed for distributed MMSE: %.2f ms.\n', td * 1000);

%performance check
if (tc < td)
    disp('Cumulative is faster than distributed!');
else 
    disp('Distributed is faster than cumulative!');
end

%sanity check
%computing relative errors
rel_errC = abs(X - x_est_c);
rel_errD = abs(X - x_est_d);

%computing mean error for x
mean_errC = mean(rel_errC);
mean_errD = mean(rel_errD);

%printing results
fprintf('\nOriginal value for x: [%.4f %.4f]\n', X);

fprintf('\nCumulative estimation for x: [%.4f %.4f]\n', x_est_c);
fprintf('Relative error on cumulative: %.2f%% and %.2f%%\n', ...
    rel_errC(1)* 100, rel_errC(2)* 100);

fprintf('\nDistributed estimation for x: [%.4f %.4f]\n', x_est_d);
fprintf('Relative error on distributed: %.2f%% and %.2f%%\n\n', ....
    rel_errD(1)* 100, rel_errD(2)* 100);

%best estimator evaluation
err_tol = 1e-04;
if (abs(mean_errC - mean_errD) < err_tol)
    disp('They have the same precision!');
elseif (mean_errC < mean_errD)
    disp('Cumulative is safer than distributed!');
else
    disp('Distributed is safer than cumulative!');
end