%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% raylee_sensitivity.m
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
% Program raylee_sensitivity is a Matlab function to compute phase and 
% group velocity sensitivity kernels of Rayleigh waves over a 
% range of frequencies.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% Nn            number of elements in solid part of model
% Nnf           number of elements in fluid part of model
% hv            vector of grid spacings in solid (meters)
% hfv           vector of grid spacings in fluid (meters)
% fks           vector of frequencies (Hz)
% modnv         vector of mode numbers (1=fundamental)
% vsv           shear velocity model, a vector (m/s)
% vpv           compressional velocity model in solid, a vector (m/s)
% rhov          density model in solid, a vector (kg/m^3)
% vpfv          compressional velocity model in fluid, a vector (m/s)
% rhofv         density model in fluid, a vector (kg/m^3)
% vflg          vector of phase or group flag (=0 for phase, =1 for group)
% pratioflag    flag indicating if P-wave velocity (=0) or 
%               Poisson's ratio (=1) is fixed
%
% Output:
% U             modeled velocities (group or phase depending on vflg) over 
%               the entire frequency range
% snsmf_vstotf  group or phase velocity sensitivity kernel (again, 
%               depending on vflg)
% snsmf_htotf   sensitivity kernel for phase  or group velocity due to an interface
%               within the layering changing its depth
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uf, snsmf_vstotf, snsmf_htotf] = ...
    raylee_sensitivity(Nn,vsv,vpv,rhov,fks,hv,modnv,vflg,Nnf,vpfv,rhofv,hfv,pratioflag)

countr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% augment the frequency vector if group kernels are needed, this triples
% the size of the frequency vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fks_orig = fks;
modnv_orig = modnv;
if (vflg(1) == 1)
    fks = [ fks_orig(1)*.999 fks_orig(1) fks_orig(1)*1.001];
    modnv = [ modnv_orig(1) modnv_orig(1) modnv_orig(1)];
else
    fks = fks_orig(1);
    modnv = modnv_orig(1);
end
for ii=2:length(fks_orig)
    
    if (vflg(ii) == 1)
        fks = [fks fks_orig(ii)*.999 fks_orig(ii) fks_orig(ii)*1.001];
        modnv = [ modnv modnv_orig(ii) modnv_orig(ii) modnv_orig(ii)];        
    else
        fks = [fks fks_orig(ii)];
        modnv = [ modnv modnv_orig(ii)];        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate phase velocities and eigenfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vp = zeros(1,length(fks));
U = zeros(1,length(fks));
xm = zeros(2*Nn,length(fks));

for f=fks

countr = countr + 1;
[kk, vpk, vgk, x] = raylee_lysmer(Nn,vsv,vpv,rhov,f,hv,modnv(countr),Nnf,vpfv,rhofv,hfv);

vp(countr) = vpk;
U(countr) = vgk;
xm(:,countr) = x;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct phase sensitivity kernels
%
% these are derivative matrices and are extremely sparse, so the 
% necessary matrix-vector multiplications are hardcoded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snsmf = zeros(Nn,length(fks));
snsmflam = zeros(Nn,length(fks));
snsmfrho = zeros(Nn,length(fks));

% for all depths (except the bottom element) and frequencies
for ii=1:(Nn-1)
    
    h = hv(ii);

    countr = 0;
for f=fks
    
    countr = countr + 1;
    
    % density sensitivity
    snsmfrho(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(2*ii-1,countr)*(h/2)*xm(2*ii-1,countr)+ ...
                            xm(2*ii+0,countr)*(h/2)*xm(2*ii+0,countr)+ ...
                            xm(2*ii+1,countr)*(h/2)*xm(2*ii+1,countr)+ ...
                            xm(2*ii+2,countr)*(h/2)*xm(2*ii+2,countr))*...
                            rhov(ii);                                                                 
    
    % lambda sensitivity
    % k^2 term
    snsmflam(ii,countr) = (1/(2*vp(countr)*U(countr)))*...
                          (xm(2*ii-1,countr)*(h/3)*xm(2*ii-1,countr)+ ...
                           xm(2*ii+1,countr)*(h/3)*xm(2*ii+1,countr)+ ...
                           xm(2*ii+1,countr)*(h/6)*xm(2*ii-1,countr)+ ...
                           xm(2*ii-1,countr)*(h/6)*xm(2*ii+1,countr))*...
                           (rhov(ii)*((vpv(ii)^2)-(2*(vsv(ii)^2))));
    %k^1 term                   
    snsmflam(ii,countr) = snsmflam(ii,countr) + ...
                          (1/(2*(2*pi*f)*U(countr)))*...
                          (xm(2*ii-1,countr)*(1/2)*xm(2*ii,countr)+ ...
                           xm(2*ii-1,countr)*(-1/2)*xm(2*ii+2,countr)+ ...
                           xm(2*ii,countr)*(1/2)*xm(2*ii-1,countr)+ ...
                           xm(2*ii,countr)*(1/2)*xm(2*ii+1,countr)+ ...
                           xm(2*ii+1,countr)*(1/2)*xm(2*ii,countr)+ ...
                           xm(2*ii+1,countr)*(-1/2)*xm(2*ii+2,countr)+ ...
                           xm(2*ii+2,countr)*(-1/2)*xm(2*ii-1,countr)+ ...
                           xm(2*ii+2,countr)*(-1/2)*xm(2*ii+1,countr))*...
                           (rhov(ii)*((vpv(ii)^2)-(2*(vsv(ii)^2))));
    % k^0 term    
    snsmflam(ii,countr) = snsmflam(ii,countr) + ...
                          (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                          (xm(2*ii+0,countr)*(1/h)*xm(2*ii+0,countr)+ ...
                           xm(2*ii+2,countr)*(1/h)*xm(2*ii+2,countr)+ ...
                           xm(2*ii+0,countr)*(-1/h)*xm(2*ii+2,countr)+ ...
                           xm(2*ii+2,countr)*(-1/h)*xm(2*ii+0,countr))*...
                           (rhov(ii)*((vpv(ii)^2)-(2*(vsv(ii)^2))));
    
    % mu sensitivity
    % k^2 term
    snsmf(ii,countr) = (1/(2*vp(countr)*U(countr)))*...
                       (xm(2*ii-1,countr)*(2*h/3)*xm(2*ii-1,countr)+ ...
                        xm(2*ii+1,countr)*(2*h/3)*xm(2*ii+1,countr)+ ...
                        xm(2*ii+1,countr)*(h/3)*xm(2*ii-1,countr)+ ...
                        xm(2*ii-1,countr)*(h/3)*xm(2*ii+1,countr)+...
                        xm(2*ii,countr)*(h/3)*xm(2*ii,countr)+ ...
                        xm(2*ii+2,countr)*(h/3)*xm(2*ii+2,countr)+ ...
                        xm(2*ii,countr)*(h/6)*xm(2*ii+2,countr)+ ...
                        xm(2*ii+2,countr)*(h/6)*xm(2*ii,countr))*...
                        (rhov(ii)*(vsv(ii)^2));
    % k^1 term                                                
    snsmf(ii,countr) = snsmf(ii,countr) + ...
                       (1/(2*(2*pi*f)*U(countr)))*...
                       (xm(2*ii-1,countr)*(-1/2)*xm(2*ii,countr)+ ...
                        xm(2*ii-1,countr)*(-1/2)*xm(2*ii+2,countr)+ ...
                        xm(2*ii,countr)*(-1/2)*xm(2*ii-1,countr)+ ...
                        xm(2*ii,countr)*(1/2)*xm(2*ii+1,countr)+ ...
                        xm(2*ii+1,countr)*(1/2)*xm(2*ii,countr)+ ...
                        xm(2*ii+1,countr)*(1/2)*xm(2*ii+2,countr)+ ...
                        xm(2*ii+2,countr)*(-1/2)*xm(2*ii-1,countr)+ ...
                        xm(2*ii+2,countr)*(1/2)*xm(2*ii+1,countr))*...
                        (rhov(ii)*(vsv(ii)^2));
    % k^0 term                                                                        
    snsmf(ii,countr) = snsmf(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(2*ii+0,countr)*(2/h)*xm(2*ii+0,countr)+ ...
                        xm(2*ii+2,countr)*(2/h)*xm(2*ii+2,countr)+ ...
                        xm(2*ii+0,countr)*(-2/h)*xm(2*ii+2,countr)+ ...
                        xm(2*ii+2,countr)*(-2/h)*xm(2*ii+0,countr)+ ...
                        xm(2*ii-1,countr)*(1/h)*xm(2*ii-1,countr)+ ...
                        xm(2*ii+1,countr)*(1/h)*xm(2*ii+1,countr)+ ...
                        xm(2*ii-1,countr)*(-1/h)*xm(2*ii+1,countr)+ ...
                        xm(2*ii+1,countr)*(-1/h)*xm(2*ii-1,countr))*...
                        (rhov(ii)*(vsv(ii)^2));  
                    
    % thickness sensitivity
    % omega^2 term
    snsmfh(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(2*ii-1,countr)*(rhov(ii)/2)*xm(2*ii-1,countr)+ ...
                            xm(2*ii+0,countr)*(rhov(ii)/2)*xm(2*ii+0,countr)+ ...
                            xm(2*ii+1,countr)*(rhov(ii)/2)*xm(2*ii+1,countr)+ ...
                            xm(2*ii+2,countr)*(rhov(ii)/2)*xm(2*ii+2,countr))*...
                            h;   
        
    % k^2 term
    pmod = (rhov(ii)*(vpv(ii)^2));
    smod = (rhov(ii)*(vsv(ii)^2)); 
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...
                        (1/(2*vp(countr)*U(countr)))*...
                       (xm(2*ii-1,countr)*(pmod/3)*xm(2*ii-1,countr)+ ...
                        xm(2*ii+1,countr)*(pmod/3)*xm(2*ii+1,countr)+ ...
                        xm(2*ii+1,countr)*(pmod/6)*xm(2*ii-1,countr)+ ...
                        xm(2*ii-1,countr)*(pmod/6)*xm(2*ii+1,countr)+...
                        xm(2*ii,countr)*(smod/3)*xm(2*ii,countr)+ ...
                        xm(2*ii+2,countr)*(smod/3)*xm(2*ii+2,countr)+ ...
                        xm(2*ii,countr)*(smod/6)*xm(2*ii+2,countr)+ ...
                        xm(2*ii+2,countr)*(smod/6)*xm(2*ii,countr))*...
                        h;    
                    
    % k^0 term                                                                        
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(2*ii+0,countr)*(pmod)*(-1/(h^2))*xm(2*ii+0,countr)+ ...
                        xm(2*ii+2,countr)*(pmod)*(-1/(h^2))*xm(2*ii+2,countr)+ ...
                        xm(2*ii+0,countr)*(-pmod)*(-1/(h^2))*xm(2*ii+2,countr)+ ...
                        xm(2*ii+2,countr)*(-pmod)*(-1/(h^2))*xm(2*ii+0,countr)+ ...
                        xm(2*ii-1,countr)*(smod)*(-1/(h^2))*xm(2*ii-1,countr)+ ...
                        xm(2*ii+1,countr)*(smod)*(-1/(h^2))*xm(2*ii+1,countr)+ ...
                        xm(2*ii-1,countr)*(-smod)*(-1/(h^2))*xm(2*ii+1,countr)+ ...
                        xm(2*ii+1,countr)*(-smod)*(-1/(h^2))*xm(2*ii-1,countr))*...
                        h;                     
        
end

end


% special case for the bottom element
ii=Nn;
h = hv(ii);
countr = 0;
for f=fks
    
    countr = countr + 1;
    
    % density sensitivity
    snsmfrho(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(2*ii-1,countr)*(h/2)*xm(2*ii-1,countr)+ ...
                            xm(2*ii+0,countr)*(h/2)*xm(2*ii+0,countr))*...
                            rhov(ii);
    
    % lambda sensitivity
    % k^2 term
    snsmflam(ii,countr) = (1/(2*vp(countr)*U(countr)))*...
                          (xm(2*ii-1,countr)*(h/3)*xm(2*ii-1,countr))*...
                          (rhov(ii)*((vpv(ii)^2)-(2*(vsv(ii)^2))));
    
    %k^1 term
    snsmflam(ii,countr) = snsmflam(ii,countr) + ...
                          (1/(2*(2*pi*f)*U(countr)))*...
                          (xm(2*ii-1,countr)*(1/2)*xm(2*ii,countr)+ ...
                           xm(2*ii,countr)*(1/2)*xm(2*ii-1,countr))*...
                           (rhov(ii)*((vpv(ii)^2)-(2*(vsv(ii)^2))));
    
    %k^0 term
    snsmflam(ii,countr) = snsmflam(ii,countr) + ...
                          (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                          (xm(2*ii+0,countr)*(1/h)*xm(2*ii+0,countr))*...
                          (rhov(ii)*((vpv(ii)^2)-(2*(vsv(ii)^2))));
        
    % mu sensitivity
    % k^2 term    
    snsmf(ii,countr) = (1/(2*vp(countr)*U(countr)))*...
                       (xm(2*ii-1,countr)*(2*h/3)*xm(2*ii-1,countr)+ ...
                        xm(2*ii,countr)*(h/3)*xm(2*ii,countr))*...
                        (rhov(ii)*(vsv(ii)^2));
   
    % k^1 term
    snsmf(ii,countr) = snsmf(ii,countr) + ...
                       (1/(2*(2*pi*f)*U(countr)))*...
                       (xm(2*ii-1,countr)*(-1/2)*xm(2*ii,countr)+ ...
                        xm(2*ii,countr)*(-1/2)*xm(2*ii-1,countr))*...
                        (rhov(ii)*(vsv(ii)^2));
         
    % k^0 term                
    snsmf(ii,countr) = snsmf(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(2*ii+0,countr)*(2/h)*xm(2*ii+0,countr)+ ...
                        xm(2*ii-1,countr)*(1/h)*xm(2*ii-1,countr))*...
                        (rhov(ii)*(vsv(ii)^2));
        
    % thickness sensitivity
    % omega^2 term
    snsmfh(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(2*ii-1,countr)*(rhov(ii)/2)*xm(2*ii-1,countr)+ ...
                            xm(2*ii+0,countr)*(rhov(ii)/2)*xm(2*ii+0,countr))*...
                            h;                    

    % k^2 term
    pmod = (rhov(ii)*(vpv(ii)^2));
    smod = (rhov(ii)*(vsv(ii)^2)); 
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...       
                        (1/(2*vp(countr)*U(countr)))*...
                       (xm(2*ii-1,countr)*(pmod/3)*xm(2*ii-1,countr)+ ...
                        xm(2*ii,countr)*(smod/3)*xm(2*ii,countr))*...
                        h;
                    
    % k^0 term                
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(2*ii+0,countr)*(pmod)*(-1/(h^2))*xm(2*ii+0,countr)+ ...
                        xm(2*ii-1,countr)*(smod)*(-1/(h^2))*xm(2*ii-1,countr))*...
                        h;                        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the shear velocity phase sensitivity kernel for frequencies of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensitivity for Vp fixed or Vp/Vs fixed
if (pratioflag == 0)
    snsmf_vs = 2*snsmf - ...
        4*transpose(transpose(snsmflam)*diag((vsv.^2)./(vpv.^2 - 2*vsv.^2)));
elseif (pratioflag == 1)
    snsmf_vs = 2*(snsmf+snsmflam);
else
end


% make a vector of the frequencies of interest
if (vflg(1) == 0)
    vfi(1) = 1;
else
    vfi(1) = 2; 
end

for ii=2:length(vflg)
    
    if (vflg(ii) == 1 & vflg(ii-1) == 0)
        
        vfi(ii) = vfi(ii-1) + 2;
        
    elseif (vflg(ii) == 1 & vflg(ii-1) == 1)
        
        vfi(ii) = vfi(ii-1) + 3;
        
    elseif (vflg(ii) == 0 & vflg(ii-1) == 0)
        
        vfi(ii) = vfi(ii-1) + 1;
        
    else
        
        vfi(ii) = vfi(ii-1) + 2;
        
    end
    
end

% compute group kernels and change relative perturbations to absolute
countr = 0;
for f=fks_orig
    
    countr = countr + 1;
    
    if (vflg(countr) == 0)
        
    Uf(countr) = vp(vfi(countr));    
    % change phase kernel for absolute perturbation in the model and data
    % instead of relative perturbation
    snsmf_vstotf(:,countr) = transpose(vp(vfi(countr))*transpose(snsmf_vs(:,vfi(countr)))*diag(vsv.^-1));

    % change phase kernel with respect to thickness to absolute perturbation
    snsmf_htotf(:,countr) = transpose(vp(vfi(countr))*transpose(snsmfh(:,vfi(countr)))*diag(hv.^-1));        
        
        
    else
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the shear velocity group sensitivity kernel,
% obtained using the method of Rodi et al. BSSA (1975)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snsmf_vstot(:,countr) = snsmf_vs(:,vfi(countr)) + ...
                        ((U(vfi(countr))/vp(vfi(countr)))*(2*pi*fks(vfi(countr)))*...
                        (snsmf_vs(:,vfi(countr)+1)-snsmf_vs(:,vfi(countr)-1))/...
                        (1*(fks(vfi(countr)+1)-fks(vfi(countr)-1))*2*pi));

% for absolute perturbations                    
snsmf_vstotf(:,countr) = transpose(vp(vfi(countr))*transpose(snsmf_vstot(:,countr))*diag(vsv.^-1));


% group velocity sensitivity kernal for changes in element thickness
snsmf_htot(:,countr) = snsmfh(:,vfi(countr)) + ...
                        ((U(vfi(countr))/vp(vfi(countr)))*(2*pi*fks(vfi(countr)))*...
                        (snsmfh(:,vfi(countr)+1)-snsmfh(:,vfi(countr)-1))/...
                        (1*(fks(vfi(countr)+1)-fks(vfi(countr)-1))*2*pi));

% for absolute perturbations                     
snsmf_htotf(:,countr) = transpose(U(vfi(countr))*transpose(snsmf_htot(:,countr))*diag(hv.^-1));

% decimate the group velocity for passback
if (isnan(U(vfi(countr)+1)) | isnan(U(vfi(countr)-1)))
    Uf(countr) = NaN;
else
    Uf(countr) = U(vfi(countr));
end
        
        
    end
    
end




