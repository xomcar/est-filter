clearvars
clc
z = tf('z');
l = tf(poly([-2, 0.5]), poly([1/3, 3, 1/4, 4]), -1);
zpk(l)
zpk(spectralFactor(l))

