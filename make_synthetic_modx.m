%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% make_synthetic_modx.m
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
% Program make_synthetic_modx is a Matlab script to model fundamental mode 
% Rayleigh wave group and phase velocities for the model described in 
% Xia et al. (1999; Geophysics). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all

% time this
tic

% grid spacings (non-uniform) in meters
h = [0.1*ones(1,200) 2*ones(1,115)];
% number of nodes
Nn = length(h);
% frequencies (Hz)
fks = [5:0.5:30];
% mode number, 1=fundamental 
modn = 1;

% no water layer on top
Nnf = 0;
vpfv = 1;
rhofv = 1;
hfv = 1;

% make a 5 layered model
layrth1 = 20; % thickness in in nodes, layrth1*h = thickness in meters
layrth2 = 23; % thickness in in nodes, layrth2*h = thickness in meters
layrth3 = 25; % thickness in in nodes, layrth3*h = thickness in meters
layrth4 = 28; % thickness in in nodes, layrth4*h = thickness in meters
layrth5 = 32; % thickness in in nodes, layrth5*h = thickness in meters

% the true model, velocity in (m/s) and density in (kg/m^3)
vplay1 = 650; vslay1 = 194; rholay1 = 1820;
vplay2 = 750; vslay2 = 270; rholay2 = 1860;
vplay3 = 1400; vslay3 = 367; rholay3 = 1910;
vplay4 = 1800; vslay4 = 485; rholay4 = 1960;
vplay5 = 2150; vslay5 = 603; rholay5 = 2020;
vplay6 = 2800; vslay6 = 740; rholay6 = 2090;

% define the material property vectors for model
vpv = [vplay1*ones(1,layrth1) vplay2*ones(1,layrth2) ...
       vplay3*ones(1,layrth3) vplay4*ones(1,layrth4) ...
       vplay5*ones(1,layrth5) ...
       vplay6*ones(1,(Nn-(layrth1+layrth2+layrth3+layrth4+layrth5)))];
vsv = [vslay1*ones(1,layrth1) vslay2*ones(1,layrth2) ...
       vslay3*ones(1,layrth3) vslay4*ones(1,layrth4) ...
       vslay5*ones(1,layrth5) ...
       vslay6*ones(1,(Nn-(layrth1+layrth2+layrth3+layrth4+layrth5)))];
rhov = [rholay1*ones(1,layrth1) rholay2*ones(1,layrth2) ...
        rholay3*ones(1,layrth3) rholay4*ones(1,layrth4) ...
        rholay5*ones(1,layrth5) ...
        rholay6*ones(1,(Nn-(layrth1+layrth2+layrth3+layrth4+layrth5)))];


% compute phase and group velocities
countr = 0;
for f=fks

[kk, vpk, vgk, ev] = ... 
    raylee_lysmer(Nn,vsv,vpv,rhov,f,h,modn,Nnf,vpfv,rhofv,hfv);

countr = countr + 1;
vp(countr) = vpk;
U(countr) = vgk;

end

% end timer
toc

% add noise to phase and group velocity 
randn('state',0);
nvctr = randn(1,length(vp));
vp = (1+0.02*nvctr).*vp;
U = (1+0.02*nvctr).*U;

% write out synthetic data
save 'modx_phase_vels.ascii' vp -ascii
save 'modx_group_vels.ascii' U -ascii
save 'modx_freqs.ascii' fks -ascii





