%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% make_initial_model_ex2.m
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
% Program make_initial_model_ex2 is a Matlab script to make an initial 
% model for use in iterative Rayleigh/Scholte wave phase or group velocity 
% inversion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output:
%
% files vp_init.txt, vs_init.txt, rho_init.txt, vpf.txt, rhof.txt, 
%       grid_values_solid.txt, grid_values_fluid.txt, input_params.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all

% low end of acceptable chi-squared window
chilo = 1;
% high end of acceptable chi-squared window
chihi = 1.5;
% maximum number of potential updates
nupdats = 16;
% model standard deviation factor
mstdfact = 2;
% smoothing scale
smscl = 1000;
% flag for if Vp/Vs ratio is fixed (1) or not (0) 
pratioflag = 1;
% number of frequency measurements
Nf = length(load('velocity_values.txt')); 

% construct a grid in solid
Nn = 240;                   % number of elements in solid
h = 250*ones(1,Nn);             % grid spacing of mesh (meters)
% construct a grid in fluid
Nnf = 10;                  % number of elements in fluid
hfv = 100*ones(1,Nnf);          % grid spacing of mesh (meters)

vpvsr = 1.7321;                 % Vp/Vs ratio
gardc = 309.6;                 % constant in Gardner relation
powr = 0.25;                  % exponent in Gardner relation

% homogeneous initial model
vslay1 = 3400; vplay1 = vslay1*vpvsr; rholay1 = gardc*(vplay1^powr);
vpv = vplay1*ones(1,Nn);
vsv = vslay1*ones(1,Nn);
rhov = rholay1*ones(1,Nn);

% make the fluid part of the model    
vplay4 = 1500; rholay4 = 1000;
% make the vectors
vpvf = vplay4*ones(1,Nnf);
rhovf = rholay4*ones(1,Nnf);    
    

% write out input parameters file
fidt = fopen('input_params.txt','w');
fprintf(fidt,'%% input parameters for Rayleigh/Scholte wave inversion\n');
fprintf(fidt,'\n');
fprintf(fidt,'%i  %% flag for fixed poisson''s ratio (0=no,1=yes)\n',pratioflag);
fprintf(fidt,'%10.5f  %% smoothness scale (m)\n',smscl);
fprintf(fidt,'%10.5f  %% a priori model standard deviation factor\n',mstdfact);
fprintf(fidt,'%i  %% maximum number of updates (iterations)\n',nupdats);
fprintf(fidt,'%i  %% number of measurements\n',Nf);
fprintf(fidt,'%i  %% number of elements in solid part of model\n',Nn);
fprintf(fidt,'%i  %% number of elements in fluid part of model\n',Nnf);
fprintf(fidt,'%10.5f  %% lower chi squared window\n',chilo);
fprintf(fidt,'%10.5f  %% higher chi squared window\n',chihi);
fclose(fidt);

% write out Vp model in single column format
if (pratioflag == 1)
    fid = fopen('vp_init.txt','w');
    for ii=1:1
        fprintf(fid,'%10.5f\n',vpvsr);
    end
    fclose(fid); 
else
    fid = fopen('vp_init.txt','w');
    for ii=1:length(vpv)
        fprintf(fid,'%10.5f\n',vpv(ii));
    end
    fclose(fid);
end

% write out Vs model in single column format
fid = fopen('vs_init.txt','w');
for ii=1:length(vsv)
    fprintf(fid,'%10.5f\n',vsv(ii));
end
fclose(fid); 

% write out density model in single column format
fid = fopen('rho_init.txt','w');
for ii=1:length(rhov)
    fprintf(fid,'%10.5f\n',rhov(ii));
end
fclose(fid); 

% write out Vp model in fluid in single column format
fid = fopen('vpf.txt','w');
for ii=1:length(vpvf)
    fprintf(fid,'%10.5f\n',vpvf(ii));
end
fclose(fid);    

% write out density model in fluid in single column format
fid = fopen('rhof.txt','w');
for ii=1:length(rhovf)
    fprintf(fid,'%10.5f\n',rhovf(ii));
end
fclose(fid); 

% write out element thicknesses in solid in single column format
fid = fopen('grid_values_solid.txt','w');
for ii=1:Nn
    fprintf(fid,'%10.5f\n',h(ii));
end
fclose(fid);

% write out element thicknesses in fluid in single column format
fid = fopen('grid_values_fluid.txt','w');
for ii=1:Nnf
    fprintf(fid,'%10.5f\n',hfv(ii));
end
fclose(fid);



