
The RAYLEE package:

A set of 11 Matlab scripts and 5 functions to forward model and invert 
a collection of Rayleigh-wave phase and/or group velocities of any order modes

Below are descriptions for the scripts, the functions, and the examples.

SCRIPTS

make_synthetic_ex1:     writes "velocity_values.txt", "velocity_values_errs.txt",
                         "frequency_values.txt", "mode_values.txt", and "vtype_values.txt"
                        calls "raylee_lysmer.m"

make_synthetic_ex2:     writes "velocity_values.txt", "velocity_values_errs.txt",
                         "frequency_values.txt", "mode_values.txt", and "vtype_values.txt"
                        calls "raylee_lysmer.m"

make_synthetic_modx:     writes "modx_phase_vels.ascii" and "modx_freqs.ascii"
                        calls "raylee_lysmer.m"

make_initial_model_ex1: writes "vp_init.txt", "vs_init.txt", "rho_init.txt",
                         "vpf.txt", "rhof.txt",
                         "grid_values_solid.txt", "grid_values_fluid.txt", 
                         and "input_params.txt"
                        no calls

make_initial_model_ex2: writes "vp_init.txt", "vs_init.txt", "rho_init.txt",
                         "vpf.txt", "rhof.txt",
                         "grid_values_solid.txt", "grid_values_fluid.txt", 
                         and "input_params.txt"
                        no calls

make_initial_model_dix: writes "velocity_values.txt", "velocity_values_errs.txt",
                         "frequency_values.txt", "mode_values.txt", "vtype_values.txt", 
                         "vp_init.txt", "vs_init.txt", "rho_init.txt",
                         "vpf.txt", "rhof.txt",
                         "grid_values_solid.txt", "grid_values_fluid.txt", 
                         and "input_params.txt"
                        no calls

raylee_invert:          reads "velocity_values.txt", "velocity_values_errs.txt",
                         "frequency_values.txt", "mode_values.txt", "vtype_values.txt", 
                         "vp_init.txt", "vs_init.txt", "rho_init.txt",
                         "vpf.txt", "rhof.txt",
                         "grid_values_solid.txt", "grid_values_fluid.txt", 
                         and "input_params.txt"
                        no write
                        calls "raylee_sensitivity.m", "raylee_lysmer.m", 
                         "linvers.m", and "check_nans.m"

plot_results_ex1:       no read/write/calls
                        to be run immediately after raylee_invert

plot_results_ex2:       no read/write/calls
                        to be run immediately after raylee_invert

plot_results_modx:      no read/write/calls
                        to be run immediately after raylee_invert

numerical_tests:        no read/write
                        calls "raylee_lysmer.m" and "raylee_sensitivity.m"


FUNCTIONS

raylee_lysmer:      INPUT
                        Nn      number of nodes in solid part of model
                        Nnf     number of nodes in fluid part of model
                        hv      vector of grid spacings for solid (meters)
                        hvf     vector of grid spacings for fluid (meters)
                        f       frequency (Hz)
                        modn    which mode (1=fundamental, 2=first overtone, etc) 
                        vsv     S-wave velocity model in solid, a vector (m/s)
                        vpv     P-wave velocity model in solid, a vector (m/s)
                        rhov    density model in solid, a vector (kg/m^3)
                        vpv     P-wave velocity model in fluid, a vector (m/s)
                        rhov    density model in fluid, a vector (kg/m^3)

                    OUTPUT
                        kk      wavenumber for the Rayleigh wave at this 
                                frequency
                        vpk     phase velocity for the Rayleigh wave at 
                                this frequency
                        vgk     group velocity for the Rayleigh wave at 
                                this frequency
                        ev      vertical and horizontal displacement 
                                eigenfunctions (mode shapes), note these 
                                are scrambled
                    DEPENDENCIES
                        "stoneley_vel.m"

raylee_sensitivity: INPUT
                        Nn           number of nodes in solid part of model
                        Nnf          number of nodes in fluid part of model
                        hv           vector of grid spacings for solid (meters)
                        hvf          vector of grid spacings for fluid (meters)
                        f            frequency (Hz)
                        modn         vector of mode numbers (1=fundamental, 2=first overtone, etc) 
                        vsv          S-wave velocity model in solid, a vector (m/s)
                        vpv          P-wave velocity model in solid, a vector (m/s)
                        rhov         density model in solid, a vector (kg/m^3)
                        vpv          P-wave velocity model in fluid, a vector (m/s)
                        rhov         density model in fluid, a vector (kg/m^3)
                        vflg         vector of phase or group flag (0=phase, 1=group)
                        pratioflag   flag indicating if P-wave velocity (=0) or 
                                     Poisson's ratio (=1) is fixed

                    OUTPUT
                        U            modeled velocities (group or phase 
                                     depending on vflg) over the entire 
                                     frequency range
                        snsmf_vstotf group or phase velocity sensitivity 
                                     kernel for Vs (again, depending on vflg)
                        snsmf_htotf  group or phase velocity sensitivity 
                                     kernel (again, depending on vflg)
                                     for an interface in the layering 
                                     changing its depth

                    DEPENDENCIES
                        "raylee_lymser.m"

linvers:            INPUT
                        U_data      velocity data to be inverted
                        U           modeled velocity data
                        snsmf_vstot the jacobian or kernel matrix
                        mcmisr      the inverse square root of the model 
                                    covariance matrix
                        dcmisr      the inverse square root of the data 
                                    covariance matrix
                        Nn          number of elements 
                        vsv         the current S-wave wave velocity model
                        vsg         the initial guess for S-wave velocity

                    OUTPUT
                        dvs         the velocity update

                    DEPENDENCIES
                        none


stoneley_vel        INPUT
                        a       Vp in solid
                        b       Vs in solid
                        c       Vp in fluid
                        f       Density in fluid
                        s       Density in solid

                    OUTPUT
                        vst     Stoneley wave velocity

                    DEPENDENCIES
                        none

check_nans          INPUT
                        U             modeled velocity data
                        U_data        velocity data to be inverted
                        fks           vector of frequencies
                        modn          vector of mode numbers
                        vflg          vector flags indicating group or phase
                        snsmf_vstotf  Vs sensitivity kernel

                    OUTPUT
                        Ur            modeled velocity data with NaNs removed
                        U_datar       velocity data to be inverted with NaNs removed
                        fksr          vector of frequencies with NaNs removed
                        fksri         vector of original frequency indicies 
                        modnr         vector of mode numbers with NaNs removed
                        vflgr         vector flags indicating group or phase with NaNs removed
                        snsmf_vstotfr Vs sensitivity kernel with NaNs removed

                    DEPENDENCIES
                        none


The codes write and read the following text files:

velocity_values.txt:      The measured Rayleigh velocities
velocity_values_errs.txt: Error bars on the measurements
frequency_values.txt:     Frequencies at which the measurements are made
mode_values.txt:          Mode number (Fundamental=1, First overtone=2, etc)
vtype_values.txt:         Vector of velocity type, either phase (0) or group (1)
vp_init.txt               Initial P-wave velocity model in solid
vs_init.txt               Initial S-wave velocity model in solid
rho_init.txt              Initial density model in solid
vpf.txt                   P-wave velocity in fluid layer (if no water layer, this is not used)
rhof.txt                  Density in fluid layer (if no water layer, this is not used)
grid_values_solid.txt     Finite element grid in solid layer
grid_values_fluid.txt     Finite element grid in fluid layer (if no water layer, this is not used)
input_params.txt          See description below

The file "input_params.txt" is automatically generated by the make_initial_model scripts. 
This file contains these quantities:

% flag for fixed poisson's ratio (0=no,1=yes)
% smoothness scale (m)
% a priori model standard deviation factor
% maximum number of updates (iterations)
% number of measurements
% number of elements in solid part of model
% number of elements in fluid part of model
% lower chi squared window
% higher chi squared window

Note that when poisson's ratio is fixed in the inversion (flag=1), the code assume initial P-wave 
velocity model contains only 1 entry which is equal to the Vp/Vs ratio. 

When there is no water layer, the number of elements in the fluid part of the model is set to 0.

The lower and higher chi-squared values define a window of chi-squared when the code can terminate
the iterative perturbational inversion. 

The smoothness scale is a regularization parameter that is typically less than the maximum sensitivity depth, 
usually one-tenth of that depth.

The a prior model standard deviation factor is usually set to 2 and this factor multiplies the average data 
standard deviation to give the model standard deviation.


EXAMPLES

Place all 11 scripts and 5 functions into the same directory. 

For the first inversion example, execute the 
following scripts in order

>> make_synthetic_ex1
>> make_initial_model_ex1
>> raylee_invert
>> plot_results_ex1

This generates Figures 2, 3, and 4 in the paper.

For the second inversion example, execute the 
following scripts in order

>> make_synthetic_ex2
>> make_initial_model_ex2
>> raylee_invert
>> plot_results_ex2

This generates Figures 5, 6, and 7 in the paper.

For the third inversion example, execute the 
following scripts in order

>> make_synthetic_modx
>> make_initial_model_dix
>> raylee_invert
>> plot_results_modx

This generates Figures 8, 9, 10, and 11 in the paper.


Matt Haney and Victor Tsai
14 July 2016




