%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% plot_results_ex2.m
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
% Program plot_results_ex2 is a Matlab script to make plots from the 
% output of program raylee_invert after running the second example. It is 
% to be run immediately after running raylee_invert.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% parameters
vpvsr = 1.7321;
% make a three layered model
layrth1 = 5; % thickness in elements, layrth1*h = thickness in meters
layrth2 = 10; % thickness in elements, layrth2*h = thickness in meters
layrth3 = 50;
% the true model
vplay1 = 4000; vslay1 = vplay1/vpvsr; 
vplay2 = 3396; vslay2 = vplay2/vpvsr; 
vplay3 = 4500; vslay3 = vplay3/vpvsr; 
vplay4 = 6000; vslay4 = vplay4/vpvsr; 
vsv_true = [vslay1*ones(1,layrth1) vslay2*ones(1,layrth2) ...
    vslay3*ones(1,layrth3) ...
       vslay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))];

% plot data comparisons   
figure
fsize = 16;
plot(fks,U_data,'bo','LineWidth',2,'MarkerSize',6); hold on
plot(fks,U,'ko','LineWidth',2,'MarkerSize',6)
plot(fksr_guess,U_guess,'ro','LineWidth',2,'MarkerSize',6)
axis([.1 .65 1500 3400])
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Velocity (m/s) '); xlabel(' Frequency (Hz) ');
title(' Data (blue), initial guess (red), and final update (black) ');

orient landscape
print(gcf,'-dpsc','Data_space_h2olyr.ps');

% plot model comparisons
figure
fsize = 16;
plot(vsv_true,0.001*hss,'b-','LineWidth',4);  
axis ij; hold on
plot(vsv_guess,0.001*hss,'r--','LineWidth',4); axis([1500 4000 0 30]); 
axis ij; hold on
plot(vsv_update((nupdat),:),0.001*hss,'k--','LineWidth',4);
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Depth (km) '); xlabel(' Shear velocity (m/s) ');
title(' True (blue) and initial models (red), and final update (black) ');

orient landscape
print(gcf,'-dpsc','Model_space_h2olyr.ps');

% sensitivity kernel of final update
figure
fsize = 14;
[qq zz] = size(snsmf_vstot);
% find 56 and 109
subplot(1,2,1)
imagesc(fks(1:56),0.001*hss,snsmf_vstot(:,1:56)./repmat(.001*h',1,56)); colormap('jet');
axis([.1 .65 0 30])
hh = colorbar('EastOutside','FontSize',fsize,...
    'FontWeight','bold','Ytick',[-0.45:0.15:0.75]);
label = sprintf(' Sensitivity (km^{-1}) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
set(gca,'Fontsize',fsize,'FontWeight','bold');
xlabel(' Frequency (Hz) '); ylabel(' Depth (km) '); 
title(' V_{S} kernel for fundamental mode ')

subplot(1,2,2)
imagesc(fks(57:109),0.001*hss,snsmf_vstot(:,57:109)./repmat(.001*h',1,53)); colormap('jet');
axis([min(fks(57:109)) max(fks(57:109)) 0 30])
hh = colorbar('EastOutside','FontSize',fsize,...
    'FontWeight','bold','Ytick',[-0.45:0.15:0.75]);
label = sprintf(' Sensitivity (km^{-1}) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
set(gca,'Fontsize',fsize,'FontWeight','bold');
xlabel(' Frequency (Hz) '); ylabel(' Depth (km) '); 
title(' V_{S} kernel for first overtone ')

orient landscape
print(gcf,'-dpsc','Kernel_h2olyr.ps');


