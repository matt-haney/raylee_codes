%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% numerical_tests.m
%
% PROGRAMMERS:
% Matt Haney and Victor Tsai
%
% Last revision date:
% 26 April 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is distributed as part of the source-code package 
%                   raylee_inversion_codes 
% that accompanies Haney and Tsai (2017). The package can be downloaded 
% from the Geophysics source-code archive at 
%                   http://software.seg.org/2017/0003/index.html
% Use of this code is subject to acceptance of the terms and conditions
% that can be found at http://software.seg.org/disclaimer.txt 
% Copyright 2017 by The Society of Exploration Geophysicists (SEG)
% Reference:
% Haney, M. M., Tsai, V. C. (2017) Perturbational and nonperturbational 
% inversion of Rayleigh-wave velocities, Geophysics, 82(3), F15-F28.
% doi: 10.1190/geo2016-0397.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program numerical_tests is a Matlab script that compares the Jacobian 
% matrix computed with the Raylee codes to Jacobian matrices shown in 
% Table 2 of Cercato (2007, GJI). It also computes the derivative of 
% phase velocity with respect to the thickness of a crustal layer using 
% perturbation theory and compares the result to the brute force method 
% of making the crust slightly thinner and slightly thicker and taking  
% finite differences.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all

% grid spacings (non-uniform) in meters
h = 250*ones(1,1600);
% number of nodes
Nn = length(h);
% frequencies (Hz)
fks = 1./[40 30 20]; 
% mode number for each frequency
% 1=fundamental, 2=first higher mode, etc
modnv = ones(1,3);

% the Novotny crust/mantle model in Cercato
vpv = zeros(1,Nn);
vsv = zeros(1,Nn);
rhov = zeros(1,Nn);
% layer 1
% the crust: 140 elements of .25 km thickness for 35 km thick
vpv(1:140) = 6000;
vsv(1:140) = 3500;
rhov(1:140) = 2700;
% layer 2
% the mante
vpv(141:end) = 8000;
vsv(141:end) = 4500;
rhov(141:end) = 3300;
% convert m, m/s, and kg/m^3 to km, km/s, and g/cc
h = h/1000;
vpv = vpv/1000;
vsv = vsv/1000;
rhov = rhov/1000;

% phase velocity flag for each period
vflg = [0 0 0];
% no fluid layer on top of this model
Nnf = 0; vpfv = 0; rhofv = 0; hfv = 0;
[U, snsmf_vstot, snsmf_h] = ...
    raylee_sensitivity(Nn,vsv,vpv,rhov,fks,h,modnv,vflg,Nnf,vpfv,rhofv,hfv,0);

% this is the Jacobian matrix shown in Table 2 of Cercato (2007, GJI)
[sum(snsmf_vstot(1:140,:))' sum(snsmf_vstot(141:end,:))']


% layer thickness derivative by perturbational method
[snsmf_h(140,:)-snsmf_h(141,:)]

% brute force method for thickness derivative:
% decrease crustal thickness by 250 m, recompute phase velocity and take 
% finite difference approximation of dc/dh
h = 250*ones(1,1600);
% layer 1
vpv(1:139) = 6000;
vsv(1:139) = 3500;
rhov(1:139) = 2700;
% layer 2
vpv(140:end) = 8000;
vsv(140:end) = 4500;
rhov(140:end) = 3300;
% convert m, m/s, and kg/m^3 to km, km/s, and g/cc
h = h/1000;
vpv = vpv/1000;
vsv = vsv/1000;
rhov = rhov/1000;
% compute phase and group velocities
countr = 0;
for f=fks
    % fundamental mode
    modn = 1;
    % no fluid layer
    Nnf = 0; vpfv = 0; rhofv = 0; hfv = 0;
    [kk, vpk, vgk, ev] = raylee_lysmer(Nn,vsv,vpv,rhov,f,h,modn,Nnf,vpfv,rhofv,hfv);
    countr = countr + 1;
    vpp2(countr) = vpk;
end
% finite difference approximation of dc/dh
(U-vpp2)/(35-34.750)


% brute force method for thickness derivative:
% increase crustal thickness by 250 m, recompute phase velocity and take 
% finite difference approximation of dc/dh
h = 250*ones(1,1600);
% layer 1
vpv(1:141) = 6000;
vsv(1:141) = 3500;
rhov(1:141) = 2700;
% layer 2
vpv(142:end) = 8000;
vsv(142:end) = 4500;
rhov(142:end) = 3300;
% convert m, m/s, and kg/m^3 to km, km/s, and g/cc
h = h/1000;
vpv = vpv/1000;
vsv = vsv/1000;
rhov = rhov/1000;
% compute phase and group velocities
countr = 0;
for f=fks
    % fundamental mode
    modn = 1;
    % no fluid layer
    Nnf = 0; vpfv = 0; rhofv = 0; hfv = 0;
    [kk, vpk, vgk, ev] = raylee_lysmer(Nn,vsv,vpv,rhov,f,h,modn,Nnf,vpfv,rhofv,hfv);
    countr = countr + 1;
    vpp2(countr) = vpk;
end
% finite difference approximation of dc/dh
(U-vpp2)/(35-35.250)
