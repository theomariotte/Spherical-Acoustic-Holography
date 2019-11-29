% construction antenne de micros

clear; clc; 
close all

% nombre de micro

ntheta = 10;
nphi = 10;
doplot = 1;
a = 30e-2;

Rm = buildSimuSphere(ntheta,nphi,a,doplot);