%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% make_initial_model_dix.m
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
% Program make_initial_model_dix is a Matlab script to make an 
% initial model for use in iterative Rayleigh/Scholte wave phase velocity 
% inversion using the Dix method (Haney and Tsai, 2015).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: Listed below at beginning of script
%
% Output:
%
% files velocity_values.txt, velocity_values_errs.txt, 
%       frequency_values.txt, mode_values.txt, vtype_values.txt
%       vp_init.txt, vs_init.txt, rho_init.txt, vpf.txt, rhof.txt, 
%       grid_values_solid.txt, grid_values_fluid.txt, input_params.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all

% time this
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% velocity data, a row vector from low to high frequency
c = load('modx_phase_vels.ascii');
% frequency raster, a row vector from low to high frequency
fr = load('modx_freqs.ascii');
% estimates of high and low error bars
c_low = c*(1-.02);
c_hi = c*(1+.02);

% Flag for which inversion to do
% 1 = Rayleigh wave phase velocity, homogeneous form
% 2 = Rayleigh wave phase velocity, power-law form w/0.25 Poisson's ratio 
% 3 = Love wave phase velocity, power-law form 
% 4 = Rayleigh wave group velocity, homogeneous form
% 5 = Rayleigh wave phase velocity, power-law form w/0.3 Poisson's ratio
% 6 = Rayleigh wave group velocity, power-law form w/0.25 Poisson's ratio 
% 7 = Love wave group velocity, power-law form
% 8 = Rayleigh wave group velocity, power-law form w/0.3 Poisson's ratio
inv_flag_ar = 1; 
% multiplicative coefficient between wavelength and sensitivity depth
xcof = 0.5;
% number of layers per wavelength at the sensitivity depth
lrho = 20; 
% range of correlation length factors to scan over
lcormultv = [10:10:1000];
% range of model variance factors to scan over
sigmscalev = [1:.1:20];
% low value of acceptable chi-squared
chiwinlo = 1;
% high value of acceptable chi-squared
chiwinhi = 1.5;
% multiple of minumum depth for model to extend to
zminfact = 1;
% multiple of maximum depth for model to extend to
zmaxfact = 4;
% Poisson's ratio (used for homogeneous model formulation and for bulding 
% the output Vp model from the Vs model)
v = 0.45;
% Nominal power law exponent (used for group velocity power law forms)
alph = 0.3;

% These parameters are needed simply to pass output to nonlinear inversion:

% Nominal density value (used for building the output density model)
rhod = 2000; %kg/m^3
% Number of fluid layers to place above model
Nnf = 0;
% Thickness of fluid layers placed above model (if Nnf=0, not subsequently 
% used)
hfval = 200; % m
% Density of fluid layers (if Nnf=0, not subsequently used)
rhovfal = 1000; % kg/m^3
% Velocity of fluid layers (if Nnf=0, not subsequently used)
vpvfal = 1500; % m/s
% Flag whether Vp/Vs ratio is to be fixed in nonlinear inversion
pratioflag = 1;
% Max number of nonlinear updates
nupdats = 15;
% Model standard deviation factor
mstdfact = 2;
% Smoothness scale, typically less than zmax
smscl = 10; % m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build layered model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% min and max depths of exponentially increasing layer thicknesses 
zmin = min(xcof*(c./fr))*zminfact;
zmax = max(xcof*(c./fr))*zmaxfact;
% maximum number of layers
Nmax = (xcof*lrho)*log(zmax/zmin);

% layer thickness assuming the last layer between zmin and zmax is smaller
% than it should be
thkso = diff([0 zmin*exp([0:1:Nmax]/(xcof*lrho)) zmax]);
% put small layers in the shallow part above zmin
thks = [ thkso(2)*ones(1,ceil(zmin/thkso(2))) thkso(2:end) ];
% the first layer is a little thinner
thks(1) = thks(1) + (zmin - thkso(2)*ceil(zmin/thkso(2)));
z = [0 cumsum(thks) inf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scan over two regularization parameters for acceptable models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size of data and model
dsz = length(fr);
msz = length(z)-1;

% Set parameter for inversion
inv_flag = inv_flag_ar;

% set some factors that are used to define kernel for homogeneous form

% f is the factor that multiplies shear wave velocity to get Rayleigh
% wave velocity
f = sqrt((8/3) - ((-16 + 56*v - 40*(v^2))/...
    (3*(2^(2/3))*(-1 + v)*(-11 + 78*v - 123*(v^2) + 56*(v^3) + ...
    3*sqrt(3)*sqrt(-5 + 36*v - 94*(v^2) + 148*(v^3) - 165*(v^4) + ...
    112*(v^5) - 32*(v^6)))^(1/3))) + ((1/(3*(-1 + v)))*(2^(2/3))*(-11 + ...
    78*v - 123*(v^2) + 56*(v^3) + 3*sqrt(3)*sqrt(-5 + 36*v - 94*(v^2) + ...
    148*(v^3) - 165*(v^4) + 112*(v^5) - 32*(v^6)))^(1/3)));
% make sure there isn't a small imaginary part
f = real(f);
% this is the denominator of the Dix expression in Haney and Tsai (2015), 
% multiplied by 8
g = -((((-2 + f^2)^3)/((1 - f^2)^(3/2))) - ...
    (((2*sqrt(2)*sqrt((-2 + f^2 + 2*v - 2*(f^2)*v)/(-1 + v))*(4 - 4*v + ...
    (f^2)*(-1 + 2*v))))/((2 - 2*v + (f^2)*(-1 + 2*v)))) + ((8*(-2 + ...
    f^2)*(2 - 2*f^2 + sqrt(2 - 2*(f^2))*sqrt((-2 + f^2 + 2*v - ...
    2*(f^2)*v)/(-1 + v))))/((-1 +f^2)*(2*sqrt(1 - f^2) + ...
    sqrt(2)*sqrt((-2 + f^2 + 2*v - 2*(f^2)*v)/(-1 + v))))));

% these are factors and exponential terms used in the homogeneous 
% formulation
fctr1 = -((f^2 - 2)^2)*(8 - 8*(f^2) + f^4)/((1 - f^2)^1.5)/g;
exp1 = 2*sqrt(1 - f^2);
fctr2 = -8*(f^2 - 2)*2*(1+sqrt(1 - f^2)*sqrt(1 + ((f^2)*((2*v - 1)/...
    (2 - 2*v)))))/sqrt(1 - f^2)/g;
exp2 = sqrt(1 - f^2)+sqrt(1 + ((f^2)*((2*v - 1)/(2 - 2*v))));
fctr3 = 4*sqrt(1 + ((f^2)*((2*v - 1)/(2 - 2*v))))*(16*(v-1) + ...
    (f^2)*(f^2 - 8)*(2*v - 1))/(2 - 2*v + (f^2)*(2*v - 1))/g;
exp3 = 2*sqrt(1 + ((f^2)*((2*v - 1)/(2 - 2*v))));

% Calculate kernel functions (Haney and Tsai, 2015) 
% f_homo = Rayleigh wave phase velocity, homogeneous form
% f_homo_U = Rayleigh wave group velocity, homogeneous form
% f_ray1 = Rayleigh wave phase velocity, power law form, Poisson's ratio
% 0.25
% f_ray2 = Rayleigh wave phase velocity, power law form, Poisson's ratio
% 0.3
% f_love = Love wave phase velocity, power law form

% angular frequency, wavenumber, and wavenumber times depth
w = 2*pi*fr; k = w./c; kz = k'*z;

% Phase velocity, homogeneous form
f_homo = fctr3*exp(-exp3*kz)+fctr2*exp(-exp2*kz)+fctr1*exp(-exp1*kz);
f_homo(:,end) = 0;

% Group velocity, homogeneous form
f_homo_U = fctr3*(1-exp3*kz).*exp(-exp3*kz)+...
           fctr2*(1-exp2*kz).*exp(-exp2*kz)+...
           fctr1*(1-exp1*kz).*exp(-exp1*kz);        
f_homo_U(:,end) = 0;        

% Poisson ratio of 0.3
f_ray2 = -130.253*exp(-1.8362*kz) + 6.88812*exp(-1.7556*kz) + ...
    271.540*exp(-1.7380*kz) - 12.1077*exp(-1.6859*kz) - ...
    183.478*exp(-1.6750*kz) + 1.70238*exp(-1.6574*kz) - ...
    142.361*exp(-1.6398*kz) + 341.126*exp(-1.6053*kz) + ...
    4.76194*exp(-1.5877*kz) - 159.030*exp(-1.5356*kz);
% Poisson ratio of 0.25
f_ray1 = -103.14*exp(-1.8586*kz) + 6.1446*exp(-1.7714*kz) + ...
    217.120*exp(-1.7555*kz) - 10.312*exp(-1.7012*kz) - ...
    160.68*exp(-1.6842*kz) + 1.2856*exp(-1.6683*kz) - ...
    115.00*exp(-1.6524*kz) + 294.66*exp(-1.6140*kz) + ...
    4.0924*exp(-1.5981*kz) - 135.49*exp(-1.5438*kz);
% Love waves
f_love = -(1+0.85^2)*exp(-2*0.85*kz);

% Construct G matrix
if inv_flag == 1
    G = diff(f_homo,1,2);
elseif inv_flag == 2
    G = diff(f_ray1,1,2);
elseif inv_flag == 3
    G = diff(f_love,1,2);
elseif inv_flag == 4
    G = diff(f_homo_U,1,2);
elseif inv_flag == 5
    G = diff(f_ray2,1,2);
elseif inv_flag == 6
    G = ((1-alph)^2)*diff(f_ray1,1,2);
elseif inv_flag == 7
    G = ((1-alph)^2)*diff(f_love,1,2);
else
    G = ((1-alph)^2)*diff(f_ray2,1,2);
end

% Depths of layer tops
z_mid = z(1:(length(z)-1));
% Defines an array of data stdevs based on the hi and low values
sigma_d_ar = ((c_hi.^2-c_low.^2)/2);
% Data covariance matrix
Cd = diag(sigma_d_ar.^2);

% Create Xia model interp (0.5c/f) and midpoint distance matrix
% find best fit power law, find slope at shallow point, and 
% do the extrapolation to the surface
% if power law gives negative velocity at surface, use linear function
% if linear function gives negative velocity at surface use shallowest 
% velocity
% repeat similarly for deepest point in model 

% Power law fit
[aap bbp] = polyfit(log(xcof*c./fr),log(c/0.88),1);
% Linear fit
[aal bbl] = polyfit(xcof*c./fr,c/0.88,1);

% extrapolation to surface
if (c(end)/.88 - exp(aap(2))*aap(1)*((xcof*c(end)/fr(end))^aap(1)) > 0)
    czro = c(end)/.88 - exp(aap(2))*aap(1)*((xcof*c(end)/fr(end))^aap(1));
    xiatrap = [c/0.88 czro];
elseif (c(end)/.88 - aal(1)*xcof*c(end)/fr(end) > 0)
    czro = c(end)/.88 - aal(1)*xcof*c(end)/fr(end);
    xiatrap = [c/0.88 czro];
else
    czro = c(end)/.88;
    xiatrap = [c/0.88 czro];
end

% extrapolation to deepest point in model 
if (c(1)/.88 - exp(aap(2))*aap(1)*((xcof*c(1)/fr(1))^(aap(1)-1))*(xcof*c(1)/fr(1)-max(z_mid)) > 0)
    cdeep = c(1)/.88 - exp(aap(2))*aap(1)*((xcof*c(1)/fr(1))^(aap(1)-1))*(xcof*c(1)/fr(1)-max(z_mid));
    xiatrap = [cdeep xiatrap];
elseif (c(1)/.88 - aal(1)*(xcof*c(1)/fr(1)-max(z_mid)) > 0)
    cdeep = c(1)/.88 - aal(1)*(xcof*c(1)/fr(1)-max(z_mid));
    xiatrap = [cdeep xiatrap];
else
    cdeep = c(1)/.88;
    xiatrap = [cdeep xiatrap];
end

% interpolate onto model 
xia_int = interp1([max(z_mid) xcof*c./fr 0],xiatrap,z_mid,'linear');

% Layer top distance matrix used in model covariance matrix
z_mid_mat = zeros(length(z_mid)); 
for i=1:length(z_mid)
    for j=1:length(z_mid)
        z_mid_mat(i,j)=abs(z_mid(i)-z_mid(j));
    end
end


% counter for acceptable models
modcountr = 0;
% counter over correlation length factors
lcntr = 0;

% scan over correlation length factors
for lcormult=lcormultv   
    
    % increment counter
    lcntr = lcntr + 1;
    
    % counter over model variance factors
    mcntr = 0;
    
% scan over model variance factors    
for sigmscale=sigmscalev  
    
    % increment counter
    mcntr = mcntr + 1;

% define model correlation length
% lcormult times the median layer thickness
l_cor = median(thks)*lcormult;

% Defines model stdev based on average data stdev
sigma_m = (mean(sigma_d_ar)*sigmscale);

% model covariance matrix
Cm = sigma_m^2*exp(-z_mid_mat/l_cor);

% Creates augmented data vector and G matrix
d_aug = [Cd^-0.5*c.^2'; Cm^-0.5*xia_int.^2'];
G_aug = [Cd^-0.5*G; Cm^-0.5];

% Does inversion using augmented d and G, 
beta_sq_r2 = (G_aug'*G_aug)\G_aug'*d_aug;

% Chi-squared of model
chi_sqd(lcntr,mcntr) = transpose(c.^2' - G*beta_sq_r2)*...
                        (Cd^-1)*(c.^2' - G*beta_sq_r2)/dsz;
                    
% Does model have negative values of squared shear velocity?
mod_neg(lcntr,mcntr) = min(beta_sq_r2>0);

% If the model is acceptable, save it in a matrix
if(chi_sqd(lcntr,mcntr) > chiwinlo & chi_sqd(lcntr,mcntr) < ... 
        chiwinhi & mod_neg(lcntr,mcntr) ~= 0)
    modcountr = modcountr + 1;
    beta_sq_r2_all(:,modcountr) = beta_sq_r2;
else
end

% End loop over model standard deviation
end

lcntr
% End loop over correlation length
end

% mean of acceptable vs models
if (modcountr > 0)
    beta_sq_r2m = mean(beta_sq_r2_all')';
elseif (min(min(chi_sqd)) > chiwinhi)
    % an error message if no models found
    error('No acceptable models found: All chi-squareds above acceptable window. Consider expanding search grid over regularization parameters or expanding acceptable chi-squared range.');
elseif (max(max(chi_sqd)) < chiwinlo)
    % an error message if no models found
    error('No acceptable models found: All chi-squareds below acceptable window. Consider expanding search grid over regularization parameters or expanding acceptable chi-squared range.');
else
    % an error message if no models found
    error('No acceptable models found: All chi-squareds above or below acceptable window. Consider densifying search grid over regularization parameters')
end




% plot Chi-squared and range of acceptable models
figure
fsize = 14;
imagesc(sigmscalev,lcormultv,chi_sqd.*mod_neg); colormap('jet');
hh = colorbar('EastOutside','FontSize',fsize,...
    'FontWeight','bold'); %,'Ytick',[-0.015:0.005:0.015]);
label = sprintf(' Chi-squared ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
hold on
contour(sigmscalev,lcormultv,chi_sqd.*mod_neg,[chiwinlo chiwinhi],'w--','LineWidth',4);
set(gca,'Fontsize',fsize,'FontWeight','bold');
title('Chi-squared misfit for Dix inversion')
xlabel('Model standard deviation factor','FontSize',fsize,'FontWeight','bold')
ylabel('Model correlation length factor','FontSize',fsize,...
    'FontWeight','bold')
orient landscape
print(gcf,'-dpsc','chi_squared.ps');

% write out input files for nonlinear inversion

Nf = length(c);
Nn = length(thks);
vout = c;
voute = c_hi-c;
fks = fr;
% Mode number, 1=fundamental (only option implemented for Dix inversion) 
modnv = ones(1,length(c));
% Type of velocity measurement, 0=phase, 1=group
if (inv_flag_ar == 1 | inv_flag_ar == 2 | inv_flag_ar == 3 | ...
        inv_flag_ar == 5)
    vtypv = zeros(1,length(c));
else
    vtypv = ones(1,length(c));
end
% Output Vs velocity model
vsv = sqrt(beta_sq_r2m(1:Nn)');
% Output Vp velocity model
vpv = sqrt((2*v-2)/(2*v-1))*vsv;
% Output density model
rhov = rhod*ones(1,Nn);
vpvf = vpvfal*ones(1,Nnf);
rhovf = rhovfal*ones(1,Nnf);
h = thks;
hfv = hfval*ones(1,Nnf);

% write out velocities in single column format
fid = fopen('velocity_values.txt','w');
for ii=1:Nf
	fprintf(fid,'%10.5f\n',vout(ii));
end
fclose(fid);

% write out velocity error bars in single column format
fid = fopen('velocity_values_errs.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',voute(ii));
end
fclose(fid);

% write out frequencies in single column format
fid = fopen('frequency_values.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',fks(ii));
end
fclose(fid);

% write out variable modnv in single column format
fid = fopen('mode_values.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',modnv(ii));
end
fclose(fid);

% write out variable vtype in single column format
fid = fopen('vtype_values.txt','w');
for ii=1:Nf
    fprintf(fid,'%10.5f\n',vtypv(ii));
end
fclose(fid);

if (pratioflag == 0)
    % write out Vp model in single column format
    fid = fopen('vp_init.txt','w');
    for ii=1:length(vpv)
        fprintf(fid,'%10.5f\n',vpv(ii));
    end
    fclose(fid);
elseif (pratioflag == 1)
    fid = fopen('vp_init.txt','w');
    fprintf(fid,'%10.5f\n',sqrt((2*v-2)/(2*v-1)));
    fclose(fid);
else
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
fprintf(fidt,'%10.5f  %% lower chi squared window\n',chiwinlo);
fprintf(fidt,'%10.5f  %% higher chi squared window\n',chiwinhi);
fclose(fidt);

% end the timer
toc


