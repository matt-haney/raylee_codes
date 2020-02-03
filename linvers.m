%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% linvers.m
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
% Program linvers is a Matlab function to do weighted damped least squares 
% inversion with lsqr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% U_data        velocity data to be inverted
% U             modeled velocity data
% snsmf_vstot   the jacobian or kernel matrix
% mcmisr        the inverse square root of the model covariance matrix
% dcmisr        the inverse square root of the data covariance matrix
% Nn            number of elements 
% vsv           the current shear wave velocity model
% vsg           the initial guess for shear wave velocity
%
% Output:
% dvs           the relative velocity update
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dvs = linvers(U_data,U,snsmf_vstot,mcmisr,dcmisr,Nn,vsv,vsg)

% transpose the jacobian or kernel matrix
G = transpose(snsmf_vstot);

% the data discrepancy
duv = (U_data' - U' + G*(vsv'-vsg'));

% scale the jacobian matrix and data
Gs = dcmisr*G;
duvs = dcmisr*duv;

% invert using lsqr
dvs = lsqr([Gs; mcmisr],...
           [duvs; zeros(Nn,1)],1e-2);


