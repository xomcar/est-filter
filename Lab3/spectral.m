function tfW = spectral(tfS)
%SPECTRALFACTOR
fprintf("+++++++++ DVD +++++++++++");
z = tf('z');

[cvZ, cvP, dK] = zpkdata(tfS, 'v');
tfW = 1;
m = 0; % grado den
n = 0; % grado num
lambdaQ = dK;

% isolo zeri stabili e costruisco numeratore
for k = 1:length(cvZ)
    if cvZ(k,1) ~= 0 && abs(cvZ(k,1)) <= 1
        tfW = tfW * (z - cvZ(k,1))
        fprintf("Aggiunto zero %f\n", cvZ(k,1));
        m = m + 1;
        lambdaQ = lambdaQ * (-cvZ(k,1))^-1;
        fprintf("Aggiornata lambda2: %f\n", lambdaQ);
    end
end

% isolo poli stabili e costruisco denominatore
for j = 1:length(cvP)
    if abs(cvP(k,1)) < 1
        tfW = tfW * (z - cvP(k,1))^-1
        fprintf("Aggiunto polo %f\n", cvP(k,1));
        n = n + 1;
        lambdaQ = lambdaQ * 1/(-cvP(k,1))^-1;
        fprintf("Aggiornata lambda2: %f\n", lambdaQ);
    end
end
% faccio radice quadrata
lambda = sqrt(lambdaQ);

tfW = z^(n-m) * lambda * tfW

end

