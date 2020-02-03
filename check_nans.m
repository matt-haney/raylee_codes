%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% check_nans.m
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
% Program check_nans is a Matlab function to look for instances when
% raylee_lysmer found that no guided mode was possible (marked NaN) and 
% write out new variables without NaNs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = ...
    check_nans(U, U_data, fks, modn, vflg, snsmf_vstot)

rcntr = 0;
for ii=1:length(fks)
    
    if (isnan(U(ii)) == 0 && isnan(U_data(ii)) == 0)   
        rcntr = rcntr + 1;
        Ur(rcntr) = U(ii);
        U_datar(rcntr) = U_data(ii);
        fksr(rcntr) = fks(ii);
        fksri(rcntr) = ii;
        modnr(rcntr) = modn(ii);
        vflgr(rcntr) = vflg(ii);
        snsmf_vstotr(:,rcntr) = snsmf_vstot(:,ii);
    else
    end
end


