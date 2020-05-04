function tfCMPSF = spectralFactor(tfc)
    [z, p, k] = zpkdata(tfc, 'v');
    %remove '0' zeroes in z and save r 
    bac = z;
    z = nonzeros(z);
    r = length(bac) - length(z);
    
    %% division of zeroes and poles
    %mstable zeros
    mstbz = [];
    %unstable zeros
    unstz = [];
    %stable poles
    stbp = [];
    %marginally unstable poles
    unstp = [];
    %divide zeroes
    for i=1:length(z)
        if abs(z(i)) <= 1
            mstbz(end+1) = z(i);
        else
            unstz(end+1) = z(i);
        end
    end
    %divide poles
    for i=1:length(p)
        if abs(p(i)) < 1
            stbp(end+1) = p(i);
        else
            unstp(end+1) = p(i);
        end
    end
    
    %% computation of lambda
    num = 1; den = 1;
    for i=1:length(unstz)
        num = num * -unstz(i);
    end
    for i=1:length(unstz)
        den = den * -unstp(i);
    end
    lambda = sqrt(k * num / den);
    
    %% output
    zeroes = zeros(1, r);
    mstbz = [mstbz, zeroes];
    ratio = tf(poly(mstbz), poly(stbp), -1);
    tfCMPSF = minreal(lambda * ratio);
end