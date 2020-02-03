%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% stoneley_vel.m
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
% Program stoneley_vel is a Matlab function that calculate the fundamental
% mode velocity of the guided wave for a model of a halfspace of water over
% a halfspace of an elastic solid. This is called a Stoneley wave since its
% velocity is less than the water velocity (i.e. it is trapped in both
% directions, up and down). In contrast, a Scholte wave occurs for a finite
% water depth when the guided wave velocity is greater than in water. The 
% Stoneley wave velocity is the solution of an eighth order polynomial.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Input parameters:
% a     Vp in solid
% b     Vs in solid
% c     Vp in fluid
% f     Density in fluid
% s     Density in solid
%
% Output parameter:
% vst   Stoneley wave velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vst = stoneley_vel(a,b,c,f,s)

% coefficient on term with velocity^-16
c16 = 256*(b^16) - (512*(b^18))/(a^2) + (256*(b^20))/(a^4);

% coefficient on term with velocity^-14
c14 = -768*(b^14) + (1280*(b^16))/(a^2) - (512*(b^18))/(a^4) - ...
    (512*(b^16))/(c^2) + (1024*(b^18))/((a^2)*(c^2)) - ...
    (512*(b^20))/((a^4)*(c^2));

% coefficient on term with velocity^-12
c12 = 832*(b^12) - (1024*(b^14))/(a^2) + (256*(b^16))/(a^4) + ...
    (256*(b^16))/(c^4) - (512*(b^18))/((a^2)*(c^4)) + ...
    (256*(b^20))/((a^4)*(c^4)) + (1536*(b^14))/(c^2) - ...
    (2560*(b^16))/((a^2)*(c^2)) + (1024*(b^18))/((a^4)*(c^2)) - ...
    (64*(b^12)*(f^2))/(s^2);

% coefficient on term with velocity^-10
c10 = -416*(b^10) + (288*(b^12))/(a^2) - (768*(b^14))/(c^4) + ...
    (1280*(b^16))/((a^2)*(c^4)) - (512*(b^18))/((a^4)*(c^4)) - ...
    (1664*(b^12))/(c^2) + (2048*(b^14))/((a^2)*(c^2)) - ...
    (512*(b^16))/((a^4)*(c^2)) + (96*(b^10)*(f^2))/(s^2) + ...
    (96*(b^12)*(f^2))/((a^2)*(s^2)) + (64*(b^12)*(f^2))/((c^2)*(s^2));

% coefficient on term with velocity^-8
c8 = 112*(b^8) - (32*(b^10))/(a^2) + (832*(b^12))/(c^4) - ...
    (1024*(b^14))/((a^2)*(c^4)) + (256*(b^16))/((a^4)*(c^4)) + ...
    (832*(b^10))/(c^2) - (576*(b^12))/((a^2)*(c^2)) - ...
    (48*(b^8)*(f^2))/(s^2) - (128*(b^10)*(f^2))/((a^2)*(s^2)) - ...
    (32*(b^12)*(f^2))/((a^4)*(s^2)) - (96*(b^10)*(f^2))/((c^2)*(s^2)) - ...
    (96*(b^12)*(f^2))/((a^2)*(c^2)*(s^2));

% coefficient on term with velocity^-6
c6 = -16*(b^6) - (416*(b^10))/(c^4) + (288*(b^12))/((a^2)*(c^4)) - ...
    (224*(b^8))/(c^2) + (64*(b^10))/((a^2)*(c^2)) + ...
    (16*(b^6)*(f^2))/(s^2) + (48*(b^8)*(f^2))/((a^2)*(s^2)) + ...
    (32*(b^10)*(f^2))/((a^4)*(s^2)) + (48*(b^8)*(f^2))/((c^2)*(s^2)) + ...
    (128*(b^10)*(f^2))/((a^2)*(c^2)*(s^2)) + ...
    (32*(b^12)*(f^2))/((a^4)*(c^2)*(s^2));

% coefficient on term with velocity^-4
c4 = (b^4) + (112*(b^8))/(c^4) - (32*(b^10))/((a^2)*(c^4)) + ...
    (32*(b^6))/(c^2) + ((b^4)*(f^4))/(s^4) - (2*(b^4)*(f^2))/(s^2) - ...
    (16*(b^6)*(f^2))/((a^2)*(s^2)) - (16*(b^6)*(f^2))/((c^2)*(s^2)) - ...
    (48*(b^8)*(f^2))/((a^2)*(c^2)*(s^2)) - ...
    (32*(b^10)*(f^2))/((a^4)*(c^2)*(s^2));

% coefficient on term with velocity^-2
c2 = -((16*(b^6))/(c^4)) - (2*(b^4))/(c^2) - ...
    (2*(b^4)*(f^4))/((a^2)*(s^4)) + (2*(b^4)*(f^2))/((a^2)*(s^2)) + ...
    (2*(b^4)*(f^2))/((c^2)*(s^2)) + (16*(b^6)*(f^2))/((a^2)*(c^2)*(s^2));

% coefficient on term with velocity^0
c0 = (b^4)/(c^4) + ((b^4)*(f^4))/((a^4)*(s^4)) - ...
    (2*(b^4)*(f^2))/((a^2)*(c^2)*(s^2));

% find roots
vsi = roots([c16 c14 c12 c10 c8 c6 c4 c2 c0]);
% take inverse square of roots
vs = sqrt(1./vsi);

% break roots into pure real and complex 
countr = 0;
counti = 0;
for ii=1:length(vs)
    if(isreal(vs(ii)))
        countr = countr + 1;
        vsre(countr) = vs(ii);
    else
        counti = counti + 1;
        vsim(counti) = vs(ii);
    end
end

% find all the real ones that are less than the liquid velocity
countr = 0;
for ii=1:length(vsre)
    if (vsre(ii) <= c)
    countr = countr + 1;
    vsre2(countr) = vsre(ii);
    else
    end
end
    
% find which root best satisfies stoneley equation 
for ii=1:length(vsre2)
    v = vsre2(ii);
    stnd(ii) = sqrt((b/v)^2 - (b/c)^2)*((1 - 2*(b/v)^2)^2 - ...
        (4*(b/v)^2)*sqrt((b/v)^2 - (b/a)^2)*sqrt((b/v)^2 - 1)) + ...
        (f/s)*sqrt((b/v)^2 - (b/a)^2);
end
[aa bb] = min(stnd);
vst = vsre2(bb);



