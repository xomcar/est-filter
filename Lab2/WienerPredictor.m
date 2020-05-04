clearvars
clc
l = tf(poly([0.5, 0,0]), poly([1/3, 3, 1/4, 4]), -1);
x = linspace(1, 10, 1000);
y = log(x);
pred(l, y, 4);


function pred(tfS, Y, k)
    %extract poles and zeroes
    [n, d] = tfdata(tfS, 'v');
    %execute polynomial division
    zeroes = roots(n);
    zeroes = [zeroes; zeros(k, 1)];
    poles = roots(d);
    poles = [poles; 0];
    [quoz, rem] = deconv(poly(zeroes), poly(poles));
    %causal part
    Caus = tf(rem, d, -1)
    %anticausal part
    AntiC = tf(quoz, 1, -1)
    %generate Wiener filter
    Wiener = minreal(tf(poly(rem), poly(d), -1))
    
end