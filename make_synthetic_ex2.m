%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% make_synthetic_ex2.m
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
% Program make_synthetic_ex2 is a Matlab script to make some synthetic 
% phase or group velocity data for Rayleigh/Scholte waves. It writes out 
% input files for eventual use by the program raylee_invert.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output:
%
% files velocity_values.txt, velocity_values_errs.txt, 
%       frequency_values.txt, mode_values.txt, vtype_values.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all

% construct a grid in solid
Nn = 240;                   % number of elements in solid
h = 250*ones(1,Nn);             % grid spacing of mesh (meters)
% construct a grid in fluid
Nnf = 10;                  % number of elements in fluid
hfv = 100*ones(1,Nnf);          % grid spacing of mesh (meters)

% construct vector of frequencies
Nf = 56;                   % number of measurements
fmin = 0.10;                    % minimum frequency (Hz)
df = 0.01;                      % frequency spacing (Hz)
fks = fmin+(df*[0:(Nf-1)]);     % vector of frequencies (Hz)
                                % these are the frequencies at which 
                                % the velocities are measured
                            
% model 2 modes, 1 fundamental and 1st higher mode 
Nf = 2*Nf;
% vector of mode numbers
modnv = [ones(1,Nf/2) 2*ones(1,Nf/2)];
% try to model over same frequencies
fks = [fks fks]; 
% what type of velocity - phase (0) or group (1)?
vtypv = zeros(1,Nf);

% medium parameters
vpvsr = 1.7321;                 % Vp/Vs ratio
gardc = 309.6;                 % constant in Gardner relation
powr = 0.25;                  % exponent in Gardner relation

% make a three layered model
layrth1 = 5; % thickness in elements, layrth1*h = thickness in meters
layrth2 = 10; % thickness in elements, layrth2*h = thickness in meters
layrth3 = 50;
% the true model in the solid
vplay1 = 4000; vslay1 = vplay1/vpvsr; rholay1 = gardc*(vplay1^powr);
vplay2 = 3396; vslay2 = vplay2/vpvsr; rholay2 = gardc*(vplay2^powr);
vplay3 = 4500; vslay3 = vplay3/vpvsr; rholay3 = gardc*(vplay3^powr);
vplay4 = 6000; vslay4 = vplay4/vpvsr; rholay4 = gardc*(vplay4^powr);
vpv = [vplay1*ones(1,layrth1) vplay2*ones(1,layrth2) ...
    vplay3*ones(1,layrth3) ...
       vplay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))];
vsv = [vslay1*ones(1,layrth1) vslay2*ones(1,layrth2) ...
    vslay3*ones(1,layrth3) ...
       vslay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))];
rhov = [rholay1*ones(1,layrth1) rholay2*ones(1,layrth2) ...
    rholay3*ones(1,layrth3) ...
        rholay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))];

% make the fluid part of the model    
vplay4 = 1500; rholay4 = 1000;
% the true model in the fluid
vpfv = vplay4*ones(1,Nnf);
rhofv = rholay4*ones(1,Nnf);    
    
% make some synthetic data
countr = 0;
for f=fks

    countr = countr + 1;
    modn = modnv(countr);
    
[kk, vpk, vgk, ev] = ... 
    raylee_lysmer(Nn,vsv,vpv,rhov,f,h,modn,Nnf,vpfv,rhofv,hfv);

if (vtypv(countr) == 0)
    vout(countr) = vpk;
else
    vout(countr) = vgk;
end

end

% add 1% Gaussian noise
% set randn to its default initial state so test is repeatable
randn('state',0);
vout = (1+0.025*randn(1,Nf)).*vout;

% take out the NaNs
countrn = 0;
Nf = length(vout)-sum(isnan(vout));
voutn = zeros(1,Nf);
fksn = voutn;
for ii=1:length(fks)
if(isnan(vout(ii)))
    %skip it
else
    countrn = countrn + 1;
    voutn(countrn) = vout(ii);
    fksn(countrn) = fks(ii);
    modnvn(countrn) = modnv(ii);
    vtypvn(countrn) = vtypv(ii);
end
end
vout = voutn;
fks = fksn;
modnv = modnvn;
vtypv = vtypvn;

% write out phase or group velocities
fid = fopen('velocity_values.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',vout(ii));
end
fclose(fid);

% write out velocity error bars in single column format
fid = fopen('velocity_values_errs.txt','w');
for ii=1:Nf
	fprintf(fid,'%10.5f\n',vout(ii)*.025);
end
fclose(fid);

% write out frequencies in single column format
fid = fopen('frequency_values.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',fks(ii));
end
fclose(fid);

% write out new variable modnv in single column format
fid = fopen('mode_values.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',modnv(ii));
end
fclose(fid);

% write out new variable vtype in single column format
fid = fopen('vtype_values.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',vtypv(ii));
end
fclose(fid);


