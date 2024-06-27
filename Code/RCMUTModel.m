function [CMUTNums] = RCMUTModel(CMUTParams)
%This function RCMUTModel() is used to simulate the behavior CMUTs with rectangular membranes, as described in our manuscript 
%"Small Signal Equivalent Circuit Model of High Performance Long Rectangular CMUT Membranes" By E.B.Dew et al. (Submitted to IEEE TUFFC)
%If this code helps you, please cite our paper.
%
%To use this program, input the details of your CMUT using the name-value parameters below.
%Each Parameter has a default value, so you only need to input the ones you wish to change.
%This function will then output a structure with your simulated results
%Note: It is usually far more convenient to create a structure of parameters, 
%then call edrc(struct), which converts the structure to a list of name-value
%parameters and then calls the RCMUTModel function with the entered parameters.
%In this documentation, we will explain each parameter, then how to interpret the output results.
%
%% Input Parameters: Device dimensions
% These dimensions are as defined in Fig. 1 of our manuscript. All parameters are in SI Units.
% a - Membrane half-width [m], Default value - 60e-6, must be numeric & positive.
% b - Membrane half-length [m], Default value - 1500e-6, must be numeric & positive.
% h - Membrane thickness [m], Default value - 5e-6, must be numeric & positive.
% d0 - Vacuum gap distance [m], Default value - 500e-9, must be numeric & positive.
% di - Insulating layer thickness [m], Default value - 360e-9, must be numeric & positive.
% p - Pressure difference between outside and inside of cavity [Pa], Default value - 101e3 - 4*1.33322e2 (1 atm - 4mTorr), must be numeric.
%
%% Input Parameters: Material properties (Membrane & Insulator)
% espr - Insulating layer relative permittivity, Default value - 3.9 (SiO2), must be numeric & positive.
% E - Membrane Young's Modulus [Pa], Default value - 148e9 (Si), must be numeric & positive.
% rho - Membrane Density [kg/m^3], Default value - 2329 (Si), must be numeric & positive.
%
%% Input Parameters: Metal electrode layer (Dimension & Properties)
% h - Electrode layer thickness [m], Default value - 250e-9, must be numeric & non-negative.
% E_metal - Electrode layer Young's Modulus [Pa], Default value - 76e9 (Au), must be numeric & positive.
% rho_metal - Electrode layer Density [kg/m^3], Default value - 19300  (Au), must be numeric & positive.
%
%% Input Parameters: CMUT Architectures with posts
% These dimensions are as defined in Fig. 11 of our manuscript and Table IV. All parameters are in SI Units.
% Device Type - Indicates which CMUT architecture to assume. Default value - "CD" ('Contiguous Dielectric' CMUT with no posts)
%   Other possible values: "EP" - (CMUT with raised 'Electrode Posts' to improve sensitivity),
%   "IIP" - (CMUT with patterned 'Isolated Isolation Posts').
%   Must be text. Defaults to CD if text not recognized and EP or IIP.
% ap - Radius of posts [m], Default value - 5e-6, must be numeric & non-negative. (EP or IIP only)
% d1 - Vacuum gap distance between membrane and posts [m], Default value - 180e-9, must be numeric & positive. d1 must be <= d0. (EP or IIP only)
% ap - Width of circular trench around posts [m], Default value - 2e-6, must be numeric & positive. (IIP only)
% n_1postRows - Number of rows with 1 central post, Default value - 58, must be integer & non-negative. (EP or IIP only)
% n_2postRows - Number of rows with 2 posts, Each post is a distance of ± xp from the center. Default value - 59, must be integer & non-negative. (EP or IIP only)
% xp - Distance between center of post and center of memrbane in 2-post rows [m], Default value - 25e-6, must be numeric & non-negative. (EP or IIP only)
%   xp must be set so that the two posts do not overlap. For EP devices this means xp > ap, for IIP devices this means xp > (ap + at).
% SetPostsAtCollapseHeight - Enabling this setting runs a static simulation with a CD device to determine the collapse point
%   Then sets d1 to equal the CD pull-in deflection.
%   This overrules the existing d1 value. Default - 0 (disabled), Enable - Set to 1. (EP or IIP only)
%
%% Input Parameters: Acoustic Medium
% Water - Toggles between air (Water = 0) or water (Water = 1) acoustic fluid medium. Default value - 0 (air medium)
%   Water = 0: Density = 1.225 kg/m^3, Speed of sound = 343 m/s.
%   Water = 1: Density = 1000 kg/m^3, Speed of sound = 1500 m/s.
%   Water = 1 also disables the secondary high-resolution frequency sweep. Must be numeric or logical.
% UseCustomFluidMedium  - Enables the parameters "Fluid Density" and "FluidSpeedOfSound" when set to 1. Default value - 0 (disabled)
% FluidDensity - Density of acoustic medium [kg/m^3], Requires UseCustomFluidMedium = 1, Default value - 920.29 (Oil), must be numeric & positive.
% FluidSpeedOfSound - Speed of sound in acoustic medium [m/s], Requires UseCustomFluidMedium = 1, Default value - 1484.1 (Oil), must be numeric & positive.
%
%% Input Parameters: Frequency Sweep
% UseDefaultFrequencySweep - Enables or disables default frequency sweep behavior. Default value - 1 (enabled).
%   In water & air, sweeps from fmin to fmax with frequency increments of fstep. 
%   In air, it will also do a finer sweep (increments fstepFine) near the predicted resonance frequency.
%   When disabled (set to 0), the frequency sweep will instead reference "CustomFrequencySweepRange". Must be numeric or logical.
% CustomFrequencySweepRange - Array of frequencies to simulate behaviour at [Hz]. Default value - [1e6:5e3:4e6].
%   Not used unless UseDefaultFrequencySweep = 0, must be numeric & positive.
% fmin - Minimum (starting) frequency of default frequency sweep [Hz], Default value - 100e3, must be numeric & positive.
% fmax - Maximum (ending) frequency of default frequency sweep [Hz], Default value - 5e6,must be numeric & positive.
% fstep - Frequency step-size (increment) of default frequency sweep [Hz], Default value - 10e3,must be numeric & positive.
% fstepFine - Frequency step-size (increment) near resonance frequency for default frequency swee [Hz], Default value - 2e3, must divide into fstep an integer number of times. (Air only)
% EnableFineSweepInWater - Enables fine sweep when Water = 1. Default value - 0 (disabled). Set to 1 to enable. Must be numeric or logical. (Water only)
% fineFsweepRangeScalingFactor - The frequency range of the fine sweep is multiplied by this number. Default value - 1. Must be positive & numeric.
%   Requires fine sweep enabled (Water = 0 or EnableFineSweepInWater = 1).
%
%% Input Parameters: Operating Points for Dynamic Model
% The dynamic model performs a frequency sweep at each operating voltage. These parameters determine the operating points:
% BiasByCollapsePercent - Toggles between setting operating points as a fraction of the collapse voltage or as a manually entered voltage
%   Default value - 1 (Operating points are determined by collapsePercent parameter), To instead use ManualBiasVoltage set to 0. Must be numeric or logical.
% collapsePercent - Array of operating points as a fraction of the collapse voltage. Requires BiasByCollapsePercent = 1.
%   Default value - [0.2 0.5 0.8 0.95] (Simulates operation @ 20%, 50%, 80% and 95% of collapse voltage.
%   All values must be greater than 0 and less than 1. (only precollapse).
% ManualBiasVoltage - Operating voltages [V]. Requires BiasByCollapsePercent = 0. Must be numeric, greater than 0, less than the collapse voltage.
%
% VAC - Sinusoidal driving signal amplitude [V] for dynamic model, must be numeric & positive.
%
%% Input Parameters: User-Experience
% SuppressPrints - Disables printing of intermediate steps. Default value - 1 (no extra printing), set to 0 to enable. Must be numeric or logical.
% PrintInputParams - Prints all function parameters. Default value - 0 (disabled), set to 1 to enable. Must be numeric or logical.
% WarningMSGBox - Gives warnings as popups instead of matlab command window. Default value - 1 (popups), set to 0 for warnings in the command window. Must be numeric or logical.
%
%% Input Parameters: Advanced
% It is usually best to leave these parameters at their default values.
% SimpleOutputs - Setting to simplify output structure, remove extra legacy outputs, Perfect for new users. Default Value - 1, set to 0 to disable. Must be numeric or logical.
% Clamped1DRadiationImpedance - Selects between Radiation impedance from 1D velocity profile and 2D square-like profile. Default - 1 (1D, recommended), 
%   set to 0 or 2 to use square-like velocity profile, which may be appropriate for lower aspect ratio devices.
%   Defaults to 1 for other values. Must be numeric.
% UsePistonRadiationImpedance - Selects between radiation impedance for clamped radiators (default) and rectangular piston (not recommended).
%   Default value - 0 (not using piston impedance). Set to 1 to use piston impedance (overrules Clamped1DRadiationImpedance). Must be numeric or logical.
% BiLayerSpringConstant - Setting that toggles whether or not the electrode stiffness is accounted for when calculating membrane stiffness. Still works when h_metal = 0. 
%   Default value - 1 (enabled), set to 0 to disable and neglect the stiffness of the electrode. Must be numeric or logical.
% ModelIsolatingTrench - Setting to model an isolating trench around bottom electrode as described in https://ieeexplore.ieee.org/document/10026669
%   Default value - 0 (disabled), set to 1 to enable, makes a Makes a trench of width 'IsolatingTrenchWidth' around the bottom electrode. Must be numeric or logical.
% IsolatingTrenchWidth - isolating width trench in [m], Default value - 5e-6, Must be numeric or logical. Requires ModelIsolatingTrench = 1.
% circularPostAreaCorrection - Scales the area of the modeled posts to be equivalent to the area of circular posts (EP or IIP only). Default value - 1 (enabled)
%   Enable: models posts as rectangles the same area as the circles, x dimensions: 2*ap, y-dimension: 2*ap*pi/4 (recommended)
%   Disable: models posts as squares with the side length 2*ap
%
%% Input Parameters: Parasitic Circuit Elements (Advanced)
% EnableSeriesResistor - Option to model a resistor in series with the voltage source in the electrical domain.
%   Default value - 0 (disabled), set to 1 to enable. Must be numeric or logical.
% Rs - Series resistance value [ohms]. Default value - 38, Unused unless EnableSeriesResistor = 1. Must be positive & numeric.
% EnableSeriesCapacitor - Option to model a capacitor in series with the voltage source in the electrical domain.
%   Default value - 0 (disabled), set to 1 to enable. Must be numeric or logical.
% Cs - Series capacitance [Farads]. Default value - 10e-9, Unused unless EnableSeriesCapacitor = 1. Must be positive & numeric.
% EnableParasiticCapacitance - Option to model a parasitic capacitor in parallel with the voltage source in the electrical domain.
%   Default value - 0 (disabled), set to 1 to enable. Must be numeric or logical.
% Cs - Parallel (parasitic) capacitance [Farads]. Default value - 20e-12, Unused unless EnableParasiticCapacitance = 1. Must be positive & numeric.
%
%% Input Parameters: Other (Advanced)
% calcHysteresis - Option to enable electrostatic modeling for decreasing voltage in addition to increasing voltage (hysteresis).
%   Enabling this setting causes the code to perform approximate calculations of the snapback voltage.
%   Default value - 0 (disabled), set to 1 to enable. Must be numeric or logical.
%   This is disabled by default because these calculations have not been validated. Enable at your own risk.
% calcFreqResponse - Option to enable calculations for dynamic model parameters, frequency response & radiation impedance.
%   Default value - 1 (enabled), set to 0 to disable. Must be numeric or logical.
%   Primary reason to disable - computation time (usually negligible)
% w0_increment - Increment in deflection calcualtion [m], Default value - 0.1e-9, must be numeric & positive & less or equal to than 1 nm.
%   small w0_increment is required for accurate electrostatic modeling, increasing this will impact accuracy.
%
%
%% Output results: Static 
% VDC - Bias voltage [V] for electrostatic analysis with increasing voltage, includes collapse points
% VINCLong (Disabled by SimpleOutputs) - Same as VDC
% VINC - Bias voltage [V] for electrostatic analysis with increasing voltage, no collapse points
% w0 - Membrane deflection [m] for electrostatic analysis with increasing voltage, includes collapse points
% w0INCLong (Disabled by SimpleOutputs) - Same as w0
% w0INC - Membrane deflection [m] for electrostatic analysis with increasing voltage, no collapse points
% C0 - CMUT capacitance [Farads] for electrostatic analysis with increasing voltage, includes collapse points
% CINCLong (Disabled by SimpleOutputs) - Same as C0
% CINC - CMUT capacitance [Farads] for electrostatic analysis with increasing voltage, no collapse points
% Vc - Collapse voltage [V]. Not applicable for post devices where membrane contacts posts via touchdown (If touchdown, disabled by SimpleOutputs)
% Vt - Touchdown voltage [V]. Voltage at which a membrane contacts the posts when not collapsed. 
%   (EP & IIP devices where membrane contacts posts via touchdown only)
%
%% Output results: Static + Hysteresis (calcHysteresis Required)
% Vsb - Snapback voltage [V].
% VDEC - Bias voltage [V] For electrostatic analysis with decreasing voltage, no collapse points unless SimpleOutputs is enabled.
% VDECLong (Disabled by SimpleOutputs) - Bias voltage [V] for electrostatic analysis with decreasing voltage, with collapse points.
% w0DEC - Membrane deflection [m] for electrostatic analysis with decreasing voltage, no collapse points unless SimpleOutputs is enabled.
% w0DECLong (Disabled by SimpleOutputs) - Membrane deflection [m] for electrostatic analysis with decreasing voltage, with collapse points.
% CDEC - CMUT capacitance [Farads] for electrostatic analysis with decreasing voltage, no collapse points unless SimpleOutputs is enabled.
% CDECLong (Disabled by SimpleOutputs) - CMUT capacitance [Farads] for electrostatic analysis with decreasing voltage, with collapse points.
%
%% Output results: Dynamic Model Parameters (calcFreqResponse Required)
% Frequency dependent quantities are in cell arrays, with an array in each cell corresponding to one operating point.
% Other parameters are arrays with one value per operating point.
% Cop - CMUT Capacitance at operating point [Farads] - element in circuit model
% wop - membrane deflection at operating point [m]
% Vop - Operating voltage
% phi - electromechanical transformer ratio in circuit model [N/V]
% Ceq - effective compliance in circuit model [m/N]
% Lm - effective mass in circuit model [kg]
%
% ZR0 - Radiation impedance in equivalent circuit model (complex & frequency dependent), [mechanical ohms], in cell array, 
% ZRADS (Disabled by SimpleOutputs) - same as ZR0
% Kax - wavenumber * a. Radiation impedance is usually calculated in terms of this quantity, in cell array, 
%
%% Output Results: Frequency Response:
% All frequency response parameters are stored in cell arrays. %Each cell contains a frequency dependent array correspond to one operating point.
% For mechanical quantities, the "out" version can be used to get the values transfered to the medium without thinking about peak vs RMS vs average.
% f - Driving frequency [Hz]
% Zin - electrical input impedance (ohms)
% G - conductance, calculated from Zin, (ohms^-1)
% Z1 (Disabled by SimpleOutputs) - electrical impedance of the mechanical circuit only (neglecting C0 and other circuit elements)
% G1 (Disabled by SimpleOutputs) - conductance of the mechanical circuit only (neglecting C0 and other circuit elements)
% fpred (Disabled by SimpleOutputs) - Predicted operating voltage for fine sweep [Hz] - requires fine sweep enabled.
% w0ac - membrane deflection (peak) [m]
% xdef (Disabled by SimpleOutputs) - same as w0ac
% v0 - membrane velocity (peak) [m/s]
% Vel (Disabled by SimpleOutputs) - same as v0
% vAvg - membrane deflection (Average) [m/s]
% vRMS - membrane deflection (RMS Average) [m/s]
% vOut - membrane deflection (RMS Average) [m/s] (same as vRMS)
% F0 - Lumped force applied to medium (peak) [N]
% Fpeak (Disabled by SimpleOutputs) - same as F0
% FAvg - Lumped force applied to medium (Average) [N]
% FRMS - Lumped force applied to medium (RMS Average) [N]
% Fout - Lumped force applied to medium (RMS Average) [N] (same as FRMS)
% P0 - Lumped pressure applied to medium (peak) [Pa]
% Ppeak (Disabled by SimpleOutputs) - same as P0
% PAvg - Lumped pressure applied to medium (Average) [Pa]
% PRMS  - Lumped pressure applied to medium (RMS Average) [Pa]
% Pout - Lumped pressure applied to medium (RMS Average) [Pa] (same as PRMS)
% P - Lumped pressure applied to medium (RMS Average) [Pa] (same as PRMS)
% PowerOut - power output (W), calculated with FRMS .* vRMS
% fr - resonance frequency (Hz)
% fa - antiresonance frequency (Hz)
% fc - "center frequency" (Hz), frequency at which the real part of Zin is maximized
% fr2 - "resonance" frequency (Hz), calculated as the frequency at which v0 is maximized, not always the same as fr.
% resFreq (Disabled by SimpleOutputs) - same as fc
% Cmeas (Disabled by SimpleOutputs) - "Measured" capacitance at each frequency, calculated from Zin [Farads]
%
% Copyright 2024 Eric Dew
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     https://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
arguments
    %Write about licensing
    % Device dimensions
    CMUTParams.a (1,1) {mustBeNumeric,mustBePositive} = 60e-6 %membrane half-width [m]
    CMUTParams.b (1,1) {mustBeNumeric,mustBePositive} = 1.5e-3 %membrane half-length [m]
    CMUTParams.h (1,1) {mustBeNumeric,mustBePositive} = 5e-6
    CMUTParams.d0 (1,1) {mustBeNumeric,mustBePositive} = 500e-9 %Vacuum cavity height [m]
    CMUTParams.di (1,1) {mustBeNumeric,mustBePositive} = 360e-9 %thickness of dielectric layer [m]
    % Material Properties
    CMUTParams.epsr (1,1) {mustBeNumeric,mustBePositive} = 3.9 %relative dielectric permittivity
    CMUTParams.E (1,1) {mustBeNumeric,mustBePositive} = 148e9 %Young's modulus, [pa]. 
    CMUTParams.rho (1,1) {mustBeNumeric,mustBePositive} = 2.329e3 %density of silicon in [kg/m^3]
    CMUTParams.p (1,1) {mustBeNumeric} = (101e3 - 4*1.33322e-4 * 1e6) %Pressure outside membrane - pressure inside membrane [Pa]. Default 1 atm - 4 mTorr.
    %Metal electrode layer
    CMUTParams.h_metal (1,1) {mustBeNumeric,mustBeNonnegative} = 250e-9;
    CMUTParams.rho_metal (1,1) {mustBeNumeric,mustBePositive} = 19.3*1e3; %kg/m^3
    CMUTParams.E_metal (1,1) {mustBeNumeric,mustBePositive} = 76e9; %Young's modulus of gold: ~ 76GPa

    %Medium Parameters:
    CMUTParams.Water (1,1) {mustBeNumericOrLogical} = 0 %determines whether air or water (immersion) is assumed. Air = 0, Water = 1.
    CMUTParams.UseCustomFluidMedium (1,1) {mustBeNumericOrLogical} = 0 % 0 (disabled) by default. Enable to use FluidDensity and FluidSpeedOfSound
    CMUTParams.FluidDensity (1,1) {mustBeNumeric,mustBePositive} = 920.29 % [kg/m^3]
    CMUTParams.FluidSpeedOfSound (1,1) {mustBeNumeric,mustBePositive} = 1484.1 % [m/s] Soybean oil @ 20C, https://iopscience.iop.org/article/10.1088/1742-6596/733/1/012040/pdf

    %Parameters related to IIP and EP devices
    CMUTParams.circularPostAreaCorrection  (1,1) {mustBeNumericOrLogical} = 1
    %determine IIP vs EP vs CD, CD default
    CMUTParams.DeviceType (1,:) {mustBeText} = 'CD'
    CMUTParams.ap (1,1) {mustBeNumeric,mustBeNonnegative} = 5e-6 %post radius
    CMUTParams.d1 (1,1) {mustBeNumeric,mustBePositive} = 180e-9 %distance to posts from undeflected membrane [m]
    CMUTParams.at (1,1) {mustBeNumeric,mustBeNonnegative} = 2e-6 %trench width [m]
    CMUTParams.n_1postRows (1,1) {mustBeInteger,mustBeNonnegative} = 58 %number of rows with 1 central post
    CMUTParams.n_2postRows (1,1) {mustBeInteger,mustBeNonnegative} = 59 %number of rows with 2 offcenter posts
    CMUTParams.xp (1,1) {mustBeNumeric,mustBePositive} = 25e-6 %location of center of the posts in the 2post rows, wrt center of the membrane, posts are at ±xp

    %% User Experience Parameters
    CMUTParams.SuppressPrints (1,1) {mustBeNumericOrLogical} = 1 %enable printing For debugging or to extract model parameters to use in a circuit simulator
    CMUTParams.PrintInputParams (1,1) {mustBeNumericOrLogical} = 0 %if you want to print all the parameters the code ran with. Also accessible with the .Params property of the output structure
    CMUTParams.WarningMSGBox (1,1) {mustBeNumericOrLogical} = 1 % (1 = enable) gives warnings in popups instead of in the MATLAB command window. (0 = disable, warnings in command window)
    CMUTParams.SimpleOutputs (1,1) {mustBeNumericOrLogical} = 1 % (1 = enable, 0 = disable) Simplifies output structure, removes extra legacy parameters, perfect for new users

    %% Analysis options:
    %These should be left at their default values unless you are doing something unusual.
    % Or if you only want the electrostatic simulation you could disable dynamic modeling
    CMUTParams.calcHysteresis (1,1) {mustBeNumericOrLogical} = 0 %calculates the snapback voltage and voltage-deflection curve for decreasing voltage as well. Method not validated.
    CMUTParams.calcFreqResponse (1,1) {mustBeNumericOrLogical} = 1 %simulates dynamic model assuming CW signals
 
    CMUTParams.UsePistonRadiationImpedance (1,1) {mustBeNumericOrLogical} = 0 %By default we use the radiation impedance of a clamped rectangular membrane. Enabling this option uses data for a piston instead.
    
    CMUTParams.Clamped1DRadiationImpedance (1,1) {mustBeNumeric} = 1 %1 - 1D velocity profile data for radiation impedance. 2 or 0 uses 2D square-like velocity profile. Otherwise - 1D, overruled by UsePistonRadiationImpedance

    CMUTParams.BiLayerSpringConstant (1,1) {mustBeNumericOrLogical} = 1 %(default: 1 - account for metal layer, 0 - neglect metal layer)

    % Bias voltages for dynamic model:
    % By default (to prevent users from accidentally trying to simulate a voltage above the collapse voltage)
    % The code simulates the response at fractions of the bias voltage (20%, 50%, 80%, 95%).
    % To change these ranges change CMUTParams.collapsePercent (must be < 1)
    % To input a bias voltage set BiasByCollapsePercent = 0 and then input
    % the range into CMUTParams.ManualBiasVoltage 
    CMUTParams.BiasByCollapsePercent (1,1) {mustBeNumericOrLogical} = 1
    CMUTParams.collapsePercent (1,:) {mustBeNumeric,  mustBeInRange(CMUTParams.collapsePercent,0,1)} = [0.2 0.5 0.8 0.95] % \% of collapse voltage, 0< <1
    %mustBeLessThanOrEqual(CMUTParams.collapsePercent, 1)
    %the ManualBiasVoltage parameter isn't used by default, have to set BBCP to 0
    CMUTParams.ManualBiasVoltage (1,:) {mustBeNumeric,mustBePositive} = [20 40 50 70 90] % in [V]

    %% Frequency sweep parameters:
    %By default this code calculates the approximate resonance frequency,
    %then does a precise sweep around that range, and a coarser sweep for frequencies outside this range to save time.
    %To sweep over some other range of frequencies, set CMUTParams.UseDefaultFrequencySweep = 0, and import the custom sweep range into CMUTParams.CustomFrequencySweepRange.
    %else, use fmin and fmax to set the maximum and minimum frequency,
    %fstep and fstepFine to set the coarse and fine sweep intervals respectively. fstep must be an integer multiple of fstepFine    
    CMUTParams.UseDefaultFrequencySweep (1,1) {mustBeNumericOrLogical} = 1 %1 uses the default frequency sweep settings, which should be sufficient for most cases, 0 uses a custom sweep
    CMUTParams.CustomFrequencySweepRange (1,:) {mustBeNumeric,mustBePositive} = [1e6:5e3:4e6] %Custom frequency range to sweep [Hz], not used unless UseDefaultFrequencySweep = 0. 
    CMUTParams.fmin (1,1) {mustBeNumeric,mustBePositive} = 100e3 %minimum frequency [Hz]
    CMUTParams.fmax (1,1) {mustBeNumeric,mustBePositive} = 5e6 %maximum frequency [Hz]
    CMUTParams.fstep (1,1) {mustBeNumeric,mustBePositive} = 10e3 %Frequency step size (coarse) [Hz]
    CMUTParams.fstepFine (1,1) {mustBeNumeric,mustBePositive} = 2e3 %Fine frequency step size near resonance frequency (air only), must divide into fstep an integer number of times [Hz]    
    CMUTParams.EnableFineSweepInWater (1,1) {mustBeNumericOrLogical} = 0 %set this to 1 to enable the fine sweep mechanics in water
    CMUTParams.fineFsweepRangeScalingFactor (1,1) {mustBeNumeric,mustBePositive} = 1 %factor to scale the size of the fine sweep range by

    %Circuit model parameters
    CMUTParams.VAC (1,1) {mustBeNumeric} = 1 %AC voltage applied to the circuit in [V]
    CMUTParams.EnableSeriesResistor (1,1) {mustBeNumericOrLogical} = 0 %Model a resistor in series with the voltage source in the electrical domain (1 - enable, 0 disable)
    CMUTParams.Rs (1,1) {mustBeNumeric,mustBePositive} = 38 %Series resistance [ohms]
    CMUTParams.EnableSeriesCapacitor (1,1) {mustBeNumericOrLogical} = 0 %Model a capacitor in series with the voltage source in the electrical domain (1 - enable, 0 disable)
    CMUTParams.Cs (1,1) {mustBeNumeric} = 10e-9 %Series capacitance [Farads]
    CMUTParams.EnableParasiticCapacitance (1,1) {mustBeNumericOrLogical} = 0 %Model a parasitic capacitor in parallel with the voltage source in the electrical domain (1 - enable, 0 disable)
    CMUTParams.Cp (1,1) {mustBeNumeric} = 20e-12 %Parallel (parasitic) capacitance [Farads]
    
    
    %Deflection increment
    CMUTParams.w0_increment (1,1) {mustBeNumeric,mustBePositive,mustBeLessThanOrEqual(CMUTParams.w0_increment,1e-9)} = 0.1e-9 %increment at which deflections are modeled [m]. I recommend against increasing this. See documentation.

    %for setting d1 to be the collapse point of a CD CMUT
    CMUTParams.SetPostsAtCollapseHeight (1,1) {mustBeNumericOrLogical} = 0
    %Model isolation trench around bottom electrode.
    CMUTParams.ModelIsolatingTrench (1,1) {mustBeNumericOrLogical} = 0 %not modeled by default.
    CMUTParams.IsolatingTrenchWidth (1,1) {mustBeNumeric,mustBeNonnegative} = 5e-6 %isolating trench in [m].

    %% undocumented
    %Optionally model a thin second layer of metal. This was created to model the thin layer of chromium under our gold electrode, but had very little impact
    %disabled by default, ie h_metal2 is set to zero. Only affects mass.
    CMUTParams.h_metal2 (1,1) {mustBeNumeric} = 0; %[m]
    CMUTParams.rho_metal2 (1,1) {mustBeNumeric} = 7.19*1e3;%kg/m^3
end %end arguments block
if (CMUTParams.IsolatingTrenchWidth >= CMUTParams.a) && CMUTParams.ModelIsolatingTrench %check for zero bottom electrode
    msg = 'Error: Isolating trench is >= the electrode width. Note that the total electrode width is 2*(a-IsolatingTrenchWidth). Please enter a valid IsolatingTrenchWidth (smaller than CMUTParams.a).';
    error(msg)
end


eps0 = 8.85418782e-12;
SuppressPrints = CMUTParams.SuppressPrints;
if CMUTParams.PrintInputParams
    disp(CMUTParams)
end
%% Process user inputs & extract parameters
% Many parameters that are used frequently are assigned to variables with shorter names
% Some contradictory user inputs are also handled here


%options
calcHyst = CMUTParams.calcHysteresis;
calcFreqResponse = CMUTParams.calcFreqResponse;
PistonZ = CMUTParams.UsePistonRadiationImpedance;

%circuit model parameters
VAC = CMUTParams.VAC; %AC Voltage
Water = CMUTParams.Water; %Water (1) or Air (0)
SeriesResistor = CMUTParams.EnableSeriesResistor; %use series resistor (1) or not (0)
Rs = CMUTParams.Rs; %series resistance value (ohms)
ParasiticCap = CMUTParams.EnableParasiticCapacitance; %use parallel parasitic capacitance (1) or not (0)
Cp  = CMUTParams.Cp; %parallel parasitic capacitance (farads)

BiasByCollapsePercent = CMUTParams.BiasByCollapsePercent; %1 = bias by %Vc, 0 = manual bias voltage
collapsePercent = CMUTParams.collapsePercent; %array of bias voltages as a % of collapse voltage
ManualBiasVoltage = CMUTParams.ManualBiasVoltage; %array of bias voltages (in V)
%the ManualBiasVoltage parameter isn't used by default, have to set BBCP to 0

%fsweep parameters
fmin = CMUTParams.fmin;
fmax = CMUTParams.fmax;
fstep = CMUTParams.fstep;
fstepP = CMUTParams.fstepFine; %precise sweeping range

if mod(fstep,fstepP) && (~Water || CMUTParams.EnableFineSweepInWater) %only complain about this if fstepP is actually used
    %fstepP must go into fstep an integer number of times
    fprintf('\nfstep should be an integer multiple of fstepFine. ')
    fstepP = fstep / ceil(fstep/fstepP);
    fprintf('Using a value of %.0f Hz for fstepFine.\n',fstepP)
end



%CMUT geometry - parameters for my membrane in [m]
a = CMUTParams.a;
b = CMUTParams.b;
h = CMUTParams.h;
d0 = CMUTParams.d0;
di = CMUTParams.di;
epsr = CMUTParams.epsr;
E = CMUTParams.E;
rho = CMUTParams.rho;
p = CMUTParams.p;
%electrostatic parameter
winc = CMUTParams.w0_increment;

if Water %medium properties
    rho_w = 1000; %kg/m^3 (water)
    c0 = 1500;
else
    rho_w = 1.225; %kg/m^3 (air)
    c0 = 343; %air
end

if CMUTParams.UseCustomFluidMedium %if enabled, use the user-inputted density and speed of sound for fluid parameters
    rho_w = CMUTParams.FluidDensity; %kg/m^3 (water)
    c0 = CMUTParams.FluidSpeedOfSound;
end

%EP/IIP Stuff:
if CMUTParams.SetPostsAtCollapseHeight %calculate "optimal" post height
    %Call this function without calling for the dynamic portion to
    %determine wPI then set d1 to wPI
    tempParams = CMUTParams;
    tempParams.calcHysteresis = 0;
    tempParams.calcFreqResponse = 0;
    tempParams.SetPostsAtCollapseHeight = 0;
    tempParams.DeviceType = 'CD';
    nvp = namedargs2cell(tempParams); %convert to a name value pair of the parameters and their values
    temp = RCMUTModel(nvp{:});

    d1 = temp.wPI;
    %disp(tempParams)
else
    d1 = CMUTParams.d1;
end

DeviceType = 0; %0 = CD, 1 = EP, 2 = IIP, Default - CD
if any(regexp(CMUTParams.DeviceType,'[cC][dD]')) || any(regexp(CMUTParams.DeviceType,'[xX][0oO]')) || any(regexp(CMUTParams.DeviceType,'[nN][oO]'))
    DeviceType = 0;

elseif any(regexp(CMUTParams.DeviceType,'[eEsS][pP]'))
    DeviceType = 1;

elseif any(regexp(CMUTParams.DeviceType,'[iI][iI][pP]')) || any(regexp(CMUTParams.DeviceType,'[sS][2345]'))
    DeviceType = 2;
end

%set ap as the post width for EP or ap+at for IIP:
if DeviceType == 2
    ap = CMUTParams.ap + CMUTParams.at;
else
    ap = CMUTParams.ap;
end

n_1postRows = CMUTParams.n_1postRows;
n_2postRows = CMUTParams.n_2postRows;
xp = CMUTParams.xp; 

%isolating trench modeling:
atr = CMUTParams.IsolatingTrenchWidth; %atr is the isolating trench width.

Sout = struct();%output structure

%% check for conflicting bottom electrode specifications:
if ((CMUTParams.a - CMUTParams.IsolatingTrenchWidth) < (xp + ap)) && CMUTParams.ModelIsolatingTrench && DeviceType == 1 && n_2postRows > 0
    %This basically checks if the trench overlaps with an EP (if both are being modeled). In that case the device doesn't really make sense.
    msg = ["Error: The specified isolating trench width overlaps with an Electrode post. This is conflicting. The minimum electrode width for this EP configuration is xp + ap.'" + ...
        "Conversely, the max electrode width allowed by this trench configuration is a - IsolatingTrenchWidth." + ...
        " To Fix this conflict, either (1), ensure that (a - IsolatingTrenchWidth) is > (xp + ap), (2) Disable either trench or EP modeling."];
    error(msg)
end

%same as above for the case where designs have 1 central EP
if ((CMUTParams.a - CMUTParams.IsolatingTrenchWidth) < (ap)) && CMUTParams.ModelIsolatingTrench && DeviceType == 1 && n_2postRows == 0 && n_1postRows > 0
    %This basically checks if the trench overlaps with an EP (if both are being modeled). In that case the device doesn't really make sense.
    msg = ["Error: The specified isolating trench width overlaps with an Electrode post. This is conflicting. The minimum electrode width for this EP configuration is ap.'" + ...
        "Conversely, the max electrode width allowed by this trench configuration is a - IsolatingTrenchWidth." + ...
        " To Fix this conflict, either (1), ensure that (a - IsolatingTrenchWidth) is > ap, (2) Disable either trench or EP modeling."];
    error(msg)
end

if ((CMUTParams.a - CMUTParams.IsolatingTrenchWidth) < (xp + ap)) && CMUTParams.ModelIsolatingTrench && DeviceType == 2 && n_2postRows > 0
    %This basically checks if the trench overlaps with an IIP (if both are being modeled). In that case the device doesn't really make sense.
    msg = ["Error: The specified isolating trench width overlaps with an isolated isolation post. This is conflicting. The minimum electrode width for this IIP configuration is xp + ap.'" + ...
        "Conversely, the max electrode width allowed by this trench configuration is a - IsolatingTrenchWidth." + ...
        " To Fix this conflict, either (1), ensure that (a - IsolatingTrenchWidth) is > (xp + ap), (2) Disable either trench or IIP modeling."];
    error(msg)
end

%same as above for the case where designs have 1 central IIP
if ((CMUTParams.a - CMUTParams.IsolatingTrenchWidth) < (ap)) && CMUTParams.ModelIsolatingTrench && DeviceType == 2 && n_2postRows == 0 && n_1postRows > 0
    %This basically checks if the trench overlaps with an EP (if both are being modeled). In that case the device doesn't really make sense.
    msg = ["Error: The specified isolating trench width overlaps with an isolated isolation post. This is conflicting. The minimum electrode width for this IIP configuration is ap.'" + ...
        "Conversely, the max electrode width allowed by this trench configuration is a - IsolatingTrenchWidth." + ...
        " To Fix this conflict, either (1), ensure that (a - IsolatingTrenchWidth) is > ap, (2) Disable either trench or IIP modeling."];
    error(msg)
end

% prevent a floating point error where xp should equal ap, but is off by a small numerical error:
if (xp < ap) && abs(xp-ap) < eps
    xp = ap;
end

if  (xp < ap) &&  DeviceType == 1 && n_2postRows > 0
    %This basically checks if the two central EP posts overlap with eachother. In that case the device doesn't really make sense.
    msg = ["Error: For the rows with two posts, the specified center location of each post (xp) is less than the width of the post (ap)." + ...
        "This means those posts are overlapping. " + ...
        "To Fix this conflict, either (1), ensure that xp > ap, (2) Set the number of rows with two posts to zero (n_2postRows = ), or (3) change the device type to 'CD'."];
    error(msg)
end

if  (xp < ap) && DeviceType == 2 && n_2postRows > 0
    %This basically checks if the two central IIP posts overlap with eachother. In that case the device doesn't really make sense.
    %Technically it would be possible to model IIPs where the center trench overlaps, but that would likely cause errors in modeling code. 
    msg = ["Error: For the rows with two posts, the specified center location of each post (xp) is less than the width of the post (ap)" + ...
        "plus the width of the trench (at). This means those posts are overlapping. " + ...
        "To Fix this conflict, either (1), ensure that xp > (ap+at), (2) Set the number of rows with two posts to zero (n_2postRows = ), or (3) change the device type to 'CD' or 'EP'."];
    error(msg)
end

if  (d1 > d0) &&  DeviceType > 0    
    % If d1 > d0, the deflection modeling, snapback modeling, etc will be wrong. Using this convention may be mathematically convenient for modeling certain
    % device types, so this may be fixed in a later update.
    msg = ["Error: The distance from the membrane to the posts (d1) is larger than the distance to the cavity bottom (d0)." + ...
        "This functionality is not currently supported. Please ensure that d1 <= d0."];
    error(msg)
end

%% Lumped modelling
%lumped parameters:
k0 = 64*b*E*h^3 / (15 * a^3); %eqn 25 in manuscript. Lumped spring constant
if CMUTParams.BiLayerSpringConstant && CMUTParams.h_metal > 0
    h_metal = CMUTParams.h_metal;
    E_ratio = E/CMUTParams.E_metal; %ratio of young's moduluses for calculation of neutral plane
    h_ratio = h/h_metal; %ratio of layer thicknesses for calculation of neutral plane
    d_np = h_metal/2 * (1 + E_ratio*h_ratio*(1+h_ratio )/(1+E_ratio*h_ratio)); %Location of the neutral plane (distance from top), Eqn 26 in manuscript
    Sout.k0Single = k0; %output the single layer k0 for debugging purposes
    k0 = b* 256/15 * 1/a^3 * (CMUTParams.E_metal * ((h_metal-d_np)^3 +d_np^3) + E * ((h_metal+h-d_np)^3 - (h_metal-d_np)^3)); %calculate k0 for bilayered beam, Eqn 28 in manuscript
    Sout.k0Bilayer = k0; %output the Bilayer k0 for debugging purposes
end
A0 = b* a* 32/15; %eqn 19 in manuscript. Effective area
 
alpha = sqrt(128/315);%alpha parameter - see eqn 31 in manuscript
beta = 8/15; %beta parameter - see eqn 32 in manuscript

m = b*4*a * (rho * h +     CMUTParams.h_metal*CMUTParams.rho_metal +  CMUTParams.h_metal2*CMUTParams.rho_metal2) ; %mass of membrane
m0 = alpha^2*m; %eqn 42 in manuscript. Effective area

f0 = sqrt(k0/(m0))/ (2*pi); %Eqn 43 - With bilayer spring constant accounted for

if ~SuppressPrints
    fprintf('\nk0: %f \n',k0);
    fprintf('A0: %.3E m^2 \n',A0);
    fprintf('f0 = sqrt(k0/(m0))/2pi =  %.3d MHz\n',f0/1e6);
end

geff = d0 + di/epsr; %effective dielectric gap
wmin = p*A0/k0;%minimum deflection
wstart = winc*ceil(wmin/winc); %second point in w0 array
if wmin>d0
    msgerr{1} = 'ERROR: Initial deflection due to atmospheric pressure is > total gap size. ';
    msgerr{2} = sprintf('Minimum deflection w0(0V) = %.2f nm. ',wmin*1e9);
    msgerr{3} = 'Ensure that the units on your dimensions and material properties are correct. If this issue persists try increasing d0.';

    if CMUTParams.WarningMSGBox
        msgbox(msgerr)
    else
        fprintf(strcat('\n',msgerr{1},'\n'))
        fprintf(strcat(msgerr{2},'\n'))
        fprintf(strcat(msgerr{3},'\n'))      
    end
    Sout.wmin = wmin;
    Sout.f0 = f0;
    Sout.k0 = k0;
    Sout.A0 = A0;
    Sout.m0 = m0;
    CMUTNums = Sout;
    return;
end

% Create an array for w0 from wmin to the max possible deflection, d0 for CD, d1 for EP, IIP
if DeviceType == 0
    w0 = [wmin wstart:winc:(d0 + winc*50)]; %include some extra points that are physically impossible to keep the derivative smooth at d0
elseif DeviceType == 1
    w0 = [wmin wstart:winc:(d1 + winc*50)]; %include some extra points that are physically impossible to keep the derivative smooth at d1
elseif DeviceType == 2
    w0 = [wmin wstart:winc:d1 (d1+winc):winc:(d0 + winc*50)]; %include some extra points that are physically impossible to keep the derivative smooth at d1 
end

if wmin>d1 && DeviceType>=1
    msgerr{1} = 'ERROR: Initial deflection due to atmospheric pressure is > d1 (available distance for membrane to deflect). ';
    msgerr{2} = sprintf('Minimum deflection w0(0V) = %.2f nm. ',wmin*1e9);
    msgerr{3} = 'Ensure that the units on your dimensions and material properties are correct. If this issue persists try increasing d1.';
    if CMUTParams.WarningMSGBox
        msgbox(msgerr)
    else
        fprintf(strcat('\n',msgerr{1},'\n'))
        fprintf(strcat(msgerr{2},'\n'))
        fprintf(strcat(msgerr{3},'\n'))
    end
    Sout.wmin = wmin;
    Sout.f0 = f0;
    Sout.k0 = k0;
    Sout.A0 = A0;
    Sout.m0 = m0;
    CMUTNums = Sout;
    return;
end


%%%%%%%%%%%%%%%%%%%%% Modelling Capacitance
gamma = w0/geff;

%Q is used to calculate the capacitance.
%Q' (QPrime) is used to calculate the capacitance derivative.
if CMUTParams.ModelIsolatingTrench
    Q = qRFun(gamma,(a-atr)/a); %Q function used for an otherwise contiguous dielectric CMUT cross seciton with a trench around it
    Qprime = qRFunPrime(gamma,(a-atr)/a);
else
    %default expression if there is no trench in the bottom electrode
    Q = qFun(gamma); %calculate q(gamma) in a separate function - eqn 13 in manuscript
    Qprime = qFunPrime(gamma); %calculate q'(gamma) in a separate function - eqn 15 in manuscript
end
switch DeviceType
    case 0 % CD CMUT:
        if CMUTParams.ModelIsolatingTrench %with isolating bottom trench:            
            l_total = 2*b - atr; %total length of the device: 2*b minus trench length
        else
            l_total = 2*b; %total length of the device: 2*b    
        end
        %total capacitance: (eqn 12 in manuscript, with l_total = 2b)
        Ct = l_total*eps0*a/geff .* Q;
        %first derivative of capacitance: (eqn 14 in manuscript, Qprime is calculated above)
        Ctprime = l_total*eps0*a/geff.^2 .* Qprime;
        %%%Second derivative of capacitance:
        %there is an analytical expression for the second derivative of Q but it's too cumbersome to be worth implementing
        %We use a numerical gradient instead. %Note: the spacing between the first and second point in w0 and gamma is different.
        Q2grad = gradient(Qprime)./gradient(gamma); %used the analytical derivative to minimize errors
        %total capacitance: second derivative
        Ct2prime = l_total*eps0*a/geff.^3 .* Q2grad;

    case 1 % EP CMUT        
        g1 = d1 + di/epsr; %effective air gap in the post region
        bstar = ap;%1/2 the y height of a single slice

        if CMUTParams.circularPostAreaCorrection 
            bstar =  bstar*pi/4; %with circular posts we make the area the same, which makes the effective thickness of each slice smaller by a factor of pi/4
        end

        %length of the 'CD' portion with no posts:
        if CMUTParams.ModelIsolatingTrench %with isolating bottom trench:
            bcd = b - n_1postRows*bstar - n_2postRows*bstar - atr/2;
            %total length is 2*bc, and there is only 1 trench along y, hence atr/2.
        else
            bcd = b - n_1postRows*bstar - n_2postRows*bstar;
        end


        gamma1 = w0/g1; %effective deflection wrt d1   
            
        %capacitance of the CD portion:
        Ccd = 2*eps0*bcd*(a)/geff .* Q;

        %capacitance of 1 slice of height 2b* with 1 EP in the center:
        C1ep = 2*eps0*bstar*(a)* (1/g1 * qRFun(gamma1,ap/a) + 1/geff *(Q - qRFun(gamma,ap/a)) );    

        %capacitance of 1 slice of height 2b* with 2 EP on the sides at position xP
        C2ep = 2*eps0*bstar*(a)*( 1/geff *qRFun(gamma,(xp-ap)/a) + 1/g1 * (qRFun(gamma1,(xp+ap)/a)  - qRFun(gamma1,(xp-ap)/a))  +  1/geff* (Q - qRFun(gamma,(xp+ap)/a) ));        

        %total capacitance  = C_cd slices + no_1postslices*C_1postslice + no2postslices *C_2postslice
        Ct = Ccd + n_1postRows*C1ep + n_2postRows*C2ep;

        %capacitance first derivative (analytical):
        % d/dw0 (Ct) = Cnp' + n1p*C1ep' + n2p*C2ep'
        CnpPrime = 2*eps0*bcd*(a)/ geff^2  * Qprime;
        C1epPrime = 2*eps0*bstar*(a) * ( 1/g1^2 * qRFunPrime(gamma1,ap/a) + ...
            1/geff^2 * ( Qprime - qRFunPrime(gamma,ap/a)  )  );
        C2epPrime =  2*eps0*bstar*(a) * ... 
            ( 1/g1^2 * ( qRFunPrime(gamma1,(ap+xp)/a) - qRFunPrime(gamma1,(xp-ap)/a)) + ...
            1/geff^2 *(Qprime - qRFunPrime(gamma,(ap+xp)/a) + qRFunPrime(gamma,(xp-ap)/a))  );
        Ctprime = CnpPrime + n_1postRows*C1epPrime + n_2postRows*C2epPrime;
      
        %numerical second derivative (no problems with numerical errors) 
        Ct2prime = gradient(Ctprime) ./ gradient(w0); %avoids numerical errors from discontinuity at first point     

        Sout.PostAreaFraction = (n_1postRows+2*n_2postRows)*(pi * ap^2) / (4*a*b); %Area fraction of EPs        
    case 2 % IIP CMUT
        bstar = ap;%1/2 the y height of a single slice

        if CMUTParams.circularPostAreaCorrection
            bstar =  bstar*pi/4; %with circular posts we make the area the same, which makes the effective thickness of each slice smaller by a factor of pi/4
        end

        %length of the 'CD' portion with no posts:
        if CMUTParams.ModelIsolatingTrench %with isolating bottom trench:
            bcd = b - n_1postRows*bstar - n_2postRows*bstar - atr/2;
            %total length is 2*bc, and there is only 1 trench along y, hence atr/2.
        else
            bcd = b - n_1postRows*bstar - n_2postRows*bstar;
        end

        %capacitance of the CD portion:
        Ccd = 2*eps0*bcd*(a)/geff .* Q;

        %capacitance of 1 slice of height 2b* with 1 IIP in the center:
        C1iip = 2*eps0*bstar*(a) /geff *(Q - qRFun(gamma,ap/a) )  ; %for an equivalent slice

        %capacitance of 1 slice of height 2b* with 2 IIP on the sides at position xP
        C2iip = 2*eps0*bstar*(a) * 1/geff * (qRFun(gamma,(xp-ap)/a) + Q - qRFun(gamma,(xp+ap)/a) ); %equivalent slice of a magical IIP with no trenches

        %total capacitance  = C_cd slices + no_1postslices*C_1postslice + no2postslices *C_2postslice
        Ct = Ccd + n_1postRows*C1iip + n_2postRows*C2iip;

        %capacitance first derivative (analytical):
        % d/dw0 (Ct) = Cnp' + n1p*C1iip' + n2p*C2iip'
        CnpPrime = 2*eps0*bcd*(a)/ geff^2  * Qprime;
        C1iipPrime = 2*eps0*bstar*(a) * ( 1/geff^2 * ( Qprime - qRFunPrime(gamma,ap/a)  )  );
        C2iipPrime =  2*eps0*bstar*(a) * ( 1/geff^2 *(Qprime - qRFunPrime(gamma,(ap+xp)/a) + qRFunPrime(gamma,(xp-ap)/a))  );
        Ctprime = CnpPrime + n_1postRows*C1iipPrime + n_2postRows*C2iipPrime;

        %numerical second derivative (no problems with numerical errors) 
        Ct2prime = gradient(Ctprime) ./ gradient(w0); %avoids numerical errors from discontinuity at first point   
        
        Sout.PostAreaFraction = (n_1postRows+2*n_2postRows)*(pi * ap^2) / (4*a*b); %Area fraction of IIPs
end

        %debug:
%         figure;plot(w0*1e9,Ct);
%         figure;plot(w0*1e9,Ctprime);
%         figure;plot(w0*1e9,Ct2prime);
%% Finding Stable operation

V = sqrt(2*  (k0*w0 - p*A0) ./ Ctprime ); %Eqn 6 in manuscript
%due to floating point errors, can get (k0*w0(1) - p*A0) = -3e-18 -> adds a complex component to V(1)
%Fix this by just setting V(1) to zero, since it should be by definition
V(1) = 0;


%collapse voltage: %occurs at max(V). Any solutions beyond the collapse point will be a lower V but unstable
%exception: in the case of touchdown, where the membrane touches the posts before the collapse voltage:
[Vc,collapseInd] = max(V); %would be better to figure out an analytical solution to this
wPI = w0(collapseInd);

%snapback voltage:
%need the index of d0 - the actual max deflection
if DeviceType == 0
    d0Index = find(w0==d0); %technically I feel like I should have a failsafe for this incase w0 doesn't contain d0
else
    d0Index = find(w0==d1); %d1 index for devices with posts
end

if isempty(d0Index) %should never get here, but sometimes weird floating point errors get here.
    wtol = winc/1e2; %this is a tolerance to get rid of floating point errors, should always be << the increment 
    if DeviceType == 0
        d0Index = find(abs(w0 -d0)<wtol);
    else
        d0Index = find(abs(w0 -d1)<wtol); %d1 index for devices with posts        
    end
    
    if isempty(d0Index) %should REALLY never get here
            disp('w0 does not contain d0');
                %could also assume it's end-50;
                %could add an extra point to w0, but then I would have to add one to V, etc. Maybe I could move this earlier in the code to prevent that.
    end   
end

if DeviceType > 0 &&wPI > d1 %In the case of EP/IIP devices that have Vc > V(w0=d1) (touchdown)
    %the previous method doesn't find Vc in the case of touchdown, it just finds a meaningless voltage.
    VcReal = sqrt(2*k0 / Ct2prime(d0Index) ); %the membrane is already in contact with the posts in this case.
    %this assumes the membrane moves no further after contacting the posts, which will definitely over-estimate the collapse voltage
    %Use the collapse voltage and collapse index in place of the touchdown voltage.

    Vc = V(d0Index);
    Vt = Vc;
    collapseInd = d0Index;

end


%handle collapse voltage; remove all unstable static solutions after the collapse voltage 
if DeviceType == 0 %CD
    %note: V(CollapseInd+1) is < Vc by definition, so just have the membrane fully collapse after a tiny bit more voltage
    deltaV = diff([V(collapseInd) V(collapseInd-1)]);
    VINC = [V(1:collapseInd) Vc+deltaV];
    w0INC = [w0(1:collapseInd) d0];
    CINC = [Ct(1:collapseInd) Ct(d0Index)];
else %EP IIP    
    if wPI > d1 %for EP/IIP devices where the collapse point is further than the posts
        VINC = V(1:d0Index);
        w0INC = [w0(1:d0Index)];
        CINC = [Ct(1:d0Index)];
        Sout.PostContactMode = 'touchdown';
    else
        deltaV = diff([V(collapseInd) V(collapseInd-1)]);
        VINC = [V(1:collapseInd) Vc+deltaV];
        w0INC = [w0(1:collapseInd) d1];
        CINC = [Ct(1:collapseInd) Ct(d0Index)];
        Sout.PostContactMode = 'snapdown';
    end
end

%Extra points to include a collapsed portion:
extrapoints = 100;

Vmax = ceil((Vc+20)/10)*10;

%adds extra points to indicate collapse
if DeviceType == 0 %CD
    w0INCLong = [w0INC d0 * ones(1,extrapoints)];
    VINCLong = [VINC linspace(Vc+2*deltaV,Vmax,extrapoints)];
    CINCLong = [CINC Ct(d0Index) *ones(1,extrapoints)];
else %EP & IIP
    w0INCLong = [w0INC d1 * ones(1,extrapoints)];
    CINCLong = [CINC Ct(d0Index) * ones(1,extrapoints)];
    if wPI > d1 %for touchdown
        ContactVs = linspace(Vc,Vmax,extrapoints+1);
        VINCLong = [VINC ContactVs(2:end)];
    else
        VINCLong = [VINC linspace(Vc+2*deltaV,Vmax,extrapoints)];
    end
end

if calcHyst
    if DeviceType == 0 %CD
        Vsb = sqrt(  2*(k0*d0 - p * A0)/ Ctprime(d0Index)); %snapback voltage
    else %EP IIP
        Vsb = sqrt(  2*(k0*d1 - p * A0)/ Ctprime(d0Index)); %snapback voltage
    end

    %decreasing voltage part of hysteresis
    %Find the w0 that the device snaps back to (w0SB)
    VDEC = V((collapseInd:-1:1)); %removed the +1
    W0SB_Ind = find((V(1:collapseInd))< Vsb); %find the indices that are less than the snapback voltage
    if isempty(W0SB_Ind) %if there is no insulator there won't be a snapback (causes an error)
        W0SB_Ind = 0; 
        if CMUTParams.calcFreqResponse
            msgerr{1} = 'Warning: Hysteresis resuslts cannot be calculated if there is no insulator. ';
            msgerr{2} = 'To avoid this error, either add an insulating layer (set d_i to some nonzero value) or set calcHyst = 0.';
            if CMUTParams.WarningMSGBox
                msgbox(msgerr)
            else
                fprintf(strcat('\n',msgerr{1},'\n'))
                fprintf(strcat(msgerr{2},'\n'))
            end
        end
    else
        W0SB_Ind = W0SB_Ind(end); %snapback deflection
    end
    
    % Creates the collapsed portion for decreasing voltage between Vc and Vsb
    snapped = ones(1,length(VDEC) - W0SB_Ind);
    if DeviceType == 0 %CD
        w0DEC = [d0*snapped w0(W0SB_Ind:-1:1)];        
    else %EP/IIP
        w0DEC = [d1*snapped w0(W0SB_Ind:-1:1)];
    end
    CDEC = [Ct(d0Index)*snapped Ct(W0SB_Ind:-1:1)]; %capacitance for decreasing voltage
    
    %adds extra points to indicate collapse
    if DeviceType == 0 %CD
        VDECLong = [linspace(Vmax,Vc,extrapoints) VDEC];
        w0DECLong = [d0 * ones(1,extrapoints) w0DEC];        
    else %EP & IIP
        if wPI > d1 %for touchdown
            VDECLong = [ContactVs(end:-1:2) VDEC];
            w0DECLong = [d1 * ones(1,extrapoints) w0DEC];
        else
            VDECLong = [linspace(Vmax,Vc,extrapoints) VDEC];
            w0DECLong = [d1 * ones(1,extrapoints) w0DEC];
        end
    end
    CDECLong = [Ct(d0Index) * ones(1,extrapoints) CDEC]; %capacitance for decreasing voltage
end

%% Load pre-calculated radiation impedance data
% Required for dynamic analysis
% Determines aspect ratio of the device and selects appropriate data based on user inputs.

if calcFreqResponse %required for any dynamic analysis. Can be omitted for speed in static case.
    %step 0: Ensure the radiation impedance data (in the data folder) is in the search path.

    CodePath = which("RCMUTModel.m"); %find the location of this file
    CodeFolder = fileparts(CodePath);
    %Then find the project folder, (up/out 1 folder), then find the data
    ProjectFolder = fullfile(CodeFolder,'..');
    DataPath = fullfile(ProjectFolder,"Data/");
    % addpath(genpath(DataPath))


    %load pre-calculated radiation impedance data:
    %step 1: calculate the closest aspect ratio - identify the file to load data from:
    q = round(b/(a),0);
    if ~SuppressPrints
        fprintf('q = %d\n',q);
    end

    %We precalculated radiation impedance for iteger aspect ratios from 1:30. If the aspect ratio is greater than 30 it warns the user here
    % However, this likely represents a very small loss in precision. See Mellow et al, for comparison of q = 30 vs q -> inf for a rectangular piston
    if q > 30
        q = 30;
        fprintf('Warning: aspect ratio q = b/a >30, using radiation impedance value for q = 30.');
    end

    alpha_Z = alpha; %Value used to convert the radiation impedance from being in terms of RMS average velocity to peak velocity
    %In the case of the 2D velocity profile v_rms/v_0 is equal to 128/315, or alpha^2. Else, the regular value of alpha is assumed.

    if PistonZ %Assume piston radiation impedance
        
        %Using piston radiation impedance data, calculated based on expressions from Mellow et al (2016).
        % RectRadZ= readmatrix(sprintf('KAvsRRvsXR_Piston_qis%d.txt',q)); %previously used text files but it's much slower to read
        load(fullfile(DataPath,'PistonZ',sprintf('KAvsRRvsXR_Piston_qis%d.mat',q)),'RectRadZ');
        %Data is formatted as follows:
        %Col 1: ka, Col 2: RR (Real part of impedance), Col 3: XR (Imaginary part of impedance)

    else %Use clamped membrane radiation impedance derived by Khorassany et al (2024).
        % select between Radiation impedance for 1D velocity profile and for 2D velocity profile using switch statement.
        switch CMUTParams.Clamped1DRadiationImpedance
            case 1
                %Load 1D radiation impedance results
                % RectRadZ= readmatrix(sprintf('KAvsRRvsXR_1D_qis%d.txt',q)); %previously used text files but it's much slower to read
                load(fullfile(DataPath,'RadiationImpedance_One_D_Velocity',sprintf('KAvsRRvsXR_1D_qis%d.mat',q)),'RectRadZ');
                %Col 1: ka, Col 2: RR, Col 3: XR
                alpha_Z = alpha;

                if q < 4
                    fprintf('Warning: Assuming 1D velocity profile with aspect ratio q = b/a < 4 may result in loss of precision. q = %i with current parameters.',q);
                end
            case {0,2}
                %Load 2D radiation impedance results
                % RectRadZ= readmatrix(sprintf('KAvsRRvsXR_2D_qis%d.txt',q)); %previously used text files but it's much slower to read
                load(fullfile(DataPath,'RadiationImpedance_Two_D_Velocity',sprintf('KAvsRRvsXR_2D_qis%d.mat',q)),'RectRadZ');
                %Col 1: ka, Col 2: RR, Col 3: XR
                alpha_Z = alpha^2;
            
            otherwise %1D
                %Load 1D radiation impedance results
                % RectRadZ= readmatrix(sprintf('KAvsRRvsXR_1D_qis%d.txt',q)); %previously used text files but it's much slower to read
                load(fullfile(DataPath,'RadiationImpedance_One_D_Velocity',sprintf('KAvsRRvsXR_1D_qis%d.mat',q)),'RectRadZ');
                %Col 1: ka, Col 2: RR, Col 3: XR
                alpha_Z = alpha;

                if q < 6
                    fprintf('Warning: Assuming 1D velocity profile with aspect ratio q = b/a < 6 may result in loss of precision. q = %i with current parameters.',q);
                end
        end % end switch statement.
    end %end non-piston radiation impedance logic

    KAX = RectRadZ(:,1);
    RadR = RectRadZ(:,2);
    RadX = RectRadZ(:,3);

    RadZ = RadR + 1i* RadX;
    RadZKa = KAX;
end %end radiation impedance loading

%% deriving circuit model paramters
if BiasByCollapsePercent
    Vop = collapsePercent*Vc;
else
    Vop = ManualBiasVoltage;
end
if ~SuppressPrints
    fprintf('\n Operating Voltages (V):\n')
    disp(Vop)
    fprintf('\n')
end

%make a structure/cell array for fsweep and frequency sweep outputs
FREQS = cell(length(Vop),1); %frequency arrays
GM_CALCs = cell(length(Vop),1); %conductance of mechanical circuit only
ZRADS = cell(length(Vop),1); % radiation impedance
Zin = cell(length(Vop),1); %input impedance
Vel = cell(length(Vop),1); %peak membrane velocity (v0)
Fpeak = cell(length(Vop),1); %Lumped force for peak circuit
xdef = cell(length(Vop),1); %w0ac 
fpred = ones(length(Vop),1); %predicted frequency for frequency sweep

Sout.Cop = zeros(length(Vop),1); %Capacitance (electrical)
Sout.Vop = zeros(length(Vop),1); %VDC
Sout.wop = zeros(length(Vop),1); %w0
Sout.phi = zeros(length(Vop),1); %gamme_E 
Sout.MechCap = zeros(length(Vop),1); %Mech Cap
Sout.MechL = zeros(length(Vop),1); %MechL

for VV = 1:length(Vop)
    if calcFreqResponse %Calculate Dynamic Model Parameters - Used to be a separate option
        if Vop(VV) > Vc
            disp('Error! Inputted V > Vc')
            break
        end

        [~,closestInd] = min(abs(VINCLong-Vop(VV))); %Want to find the closest point that we've calculated data for already
        Vop_calculated = VINCLong(closestInd); %if I use V, I can get the unstable solution
        wop = w0INCLong(closestInd);

        Cop = Ct(closestInd); %parallel capacitance

        if ParasiticCap
            Cop = Cop + Cp;
        end
       
        EM_ratio = Vop(VV)*Ctprime(closestInd); %phi
        keffop = k0 - 1/2 * Vop(VV)^2 * Ct2prime(closestInd);
        MechCap = 1/keffop;
        MechInductor = m0;     
        
        if ~SuppressPrints
            %display predicted results from simple \sqrt(keff/m) analysis
            fprintf('\nAt an operating voltage of %.2f V (%.2f%% of Vc):\n',Vop_calculated, Vop_calculated/Vc *100);
            fprintf('Electrical Domain Capacitance Value: %.2d F\n', Cop);
            fprintf('Electromechanical Turns Ratio \x2CAA: %.2d C/m\n', EM_ratio); %phi is in coulombs /m
            fprintf('Effective Spring Constant: %.2d N/m\n', keffop);
            fprintf('Mechanical Domain Capacitance Value: %.2d m/N\n', MechCap);
            fprintf('Mechanical Domain Inductor Value: %.2d kg\n', MechInductor);
            fprintf('Predicted Mechanical Resonance: %.3f MHz\n', 1/(2*pi *sqrt(MechCap*MechInductor)) /1e6);
        end
        Sout.Cop(VV) = Cop;
        Sout.Vop(VV) = Vop(VV);
        Sout.wop(VV) = wop;
        Sout.phi(VV) = EM_ratio;
        Sout.MechCap(VV) = MechCap;
        Sout.MechL(VV) = MechInductor; %changed from m0
    end %end calculate dynamic model parameters

    if calcFreqResponse
        fpred(VV) = 1/(2*pi *sqrt(MechCap*MechInductor)); %predicted resonance frequency
        
        if fpred(VV) > fmax && ~Water
            fprintf('\nPredicted resonance frequency %.3f MHz > fmax. Consider increasing fmax to obtain more accurate results.\n',fpred(VV)/1e6)
        end

        if CMUTParams.UseDefaultFrequencySweep
            %%%calculate frequency sweep range.
            if ~Water || CMUTParams.EnableFineSweepInWater %in air we do a fine sweep to get higher frequency resolution
                %determine the range of the precise sweep
                %precise range is a fraction of fpred
                if Water %Water adds complexity here. The change depends on the aspect ratio and radiation impedance                    
                    PreciseFreqFraction = 0.3; %use a broader range because it's harder to guess correctly
                    %approximate the change in center freq caused by immersion media
                    if q > 20
                        freq_change = 0.7;
                    elseif q> 10
                        freq_change = 0.75;
                    else
                        freq_change = 0.9;
                    end
                    if alpha_Z == alpha
                        freq_change = freq_change^2; %1D radiation impedance has greater effect on frequency:
                    end
                    fpred(VV) = fpred(VV)*freq_change;
                else
                    PreciseFreqFraction = 0.15;
                end
                PreciseFreqFraction = PreciseFreqFraction*CMUTParams.fineFsweepRangeScalingFactor; %make the fine range larger by the input scaling amount

                fPreciseRangeLow = floor((1-PreciseFreqFraction)*fpred(VV) /fstep)*fstep;
                fPreciseRangeHigh = ceil((1+PreciseFreqFraction)*fpred(VV) /fstep)*fstep;
                
                %ensure that the scaling amount doesn't take us outside the min and max frequency
                if fPreciseRangeLow < (fmin+fstep)
                    fPreciseRangeLow = fmin+fstep;
                end
                if fPreciseRangeHigh > (fmax-fstep)
                    fPreciseRangeHigh = fmax-fstep;
                end
                %a precise range including fpred from fPreciseRangeLow to fPreciseRangeHigh
                %coarse range above that

                if fpred(VV) < fmax
                    fsweep = [fmin:fstep:fPreciseRangeLow (fPreciseRangeLow+fstepP):fstepP:floor(fstepP*floor(fpred(VV)/fstepP)) fpred(VV) (fstepP*ceil(fpred(VV)/fstepP)):fstepP:fPreciseRangeHigh (fPreciseRangeHigh + fstep):fstep:fmax];
                else
                    fsweep = [fmin:fstep:fPreciseRangeLow (fPreciseRangeLow+fstepP):fstepP:fmax];
                end

            else %in immersion media, it is usually broadband enough that the fine sweep is unneccessary
                fsweep = fmin:fstep:fmax;
            end
           
        else %full custom sweep
            %disadvantage: Default settings pick the range close to the predicted resonance frequency which is more efficient
            fsweep = CMUTParams.CustomFrequencySweepRange;
        end

        kax = 2*pi*fsweep/c0 * a; %ka calculation

        %Interpolate radiation impedance values from pre-loaded data
        A = 4*a*b; %Total Area of the membrane:
        Rload = A*rho_w * c0 * interp1(RadZKa,real(RadZ),kax,'spline');
        Xload = A*rho_w * c0 * interp1(RadZKa,imag(RadZ),kax,'spline');

        Zrad = alpha_Z^2 * (Rload + 1i*Xload); %Radiation impedance in terms of peak velocity

        %some variables are renamed for the circuit analysis
        phi(VV) = EM_ratio;
        Ceq(VV) = MechCap;
        Lm(VV) = MechInductor;

        f_points = length(fsweep);
        G_calc = zeros(f_points,1);%preallocate
        Z_calc = zeros(f_points,1);%preallocate
        GM_calc = zeros(f_points,1);%preallocate %test of just the mechanical domain
        vel_calc = zeros(f_points,1);%preallocate %test of just the mechanical domain
        Fpeak_calc = zeros(f_points,1);%preallocate %test of just the mechanical domain
        x_calc = zeros(f_points,1);%preallocate %test of just the mechanical domain
        Z_1_calc = zeros(f_points,1);
        C_meas = zeros(f_points,1);
        

        for ii = 1:length(fsweep) %circuuit model calculations
            f=fsweep(ii);       

            %calculate z
            if SeriesResistor || CMUTParams.EnableSeriesCapacitor %analysis with a resistor in the electrical domain in series with the voltage source
                Z_c0d = -1i*(1/(2*pi*f * Cop));

                %impedance of the mechanical circuit in the electrical domain
                Z_1 = 1/phi(VV)^2 * (-1i * 1/(2*pi*f*Ceq(VV)) + 1i*(2*pi*f)*Lm(VV) + Zrad(ii));               
                
                if SeriesResistor && CMUTParams.EnableSeriesCapacitor %if there is both Rs and Cs
                    Zcs = -1i*(1/(2*pi*f * CMUTParams.Cs));
                    Z_calc(ii) = Rs + Zcs + Z_c0d*Z_1 / (Z_c0d+Z_1);
                elseif SeriesResistor && ~CMUTParams.EnableSeriesCapacitor  %if Rs and not Cs
                    Z_calc(ii) = Rs + Z_c0d*Z_1 / (Z_c0d+Z_1);
                else %if Cs but not Rs
                    Zcs = -1i*(1/(2*pi*f * CMUTParams.Cs));
                    Z_calc(ii) = Zcs + Z_c0d*Z_1 / (Z_c0d+Z_1);
                end                
                G_calc(ii) = real(1/Z_calc(ii));
                GM_calc(ii)= real(1/Z_1);
                I_total = VAC ./ Z_calc(ii); %total current entering the electrical port
                I_Z1 = I_total * Z_c0d ./ (Z_1 + Z_c0d); %current divider. current entering the electrical transformer
                vel_calc(ii) = I_Z1 /phi(VV);
                Fpeak_calc(ii) = -vel_calc(ii) * Zrad(ii); %added negative sign because that's how it's defined in koymen model
                Z_1_calc(ii) = Z_1;
            else
                Z_c0d = -1i*(1/(2*pi*f * Cop));
                %impedance of the mechanical circuit in the electrical domain
                Z_1 = 1/phi(VV)^2 * (-1i * 1/(2*pi*f*Ceq(VV)) + 1i*(2*pi*f)*Lm(VV) + Zrad(ii));
                Z_calc(ii) = Z_c0d*Z_1 / (Z_c0d+Z_1);
                G_calc(ii) = real(1/Z_calc(ii));
                GM_calc(ii)= real(1/Z_1);
                vel_calc(ii) = VAC / Z_1 /phi(VV);
                Fpeak_calc(ii) = -vel_calc(ii) * Zrad(ii);
                Z_1_calc(ii) = Z_1;
            end
            x_calc(ii) = 1/(1i*2*pi*f)*vel_calc(ii); 
            C_meas(ii) = -real(1./(2*pi*f*imag(Z_calc(ii))) ); %measured/total capacitance
            %C_meas needs a negative sign because of how the imaginary 1is cancel
        end %end fsweep for loop
        GM_CALCs{VV} = GM_calc;
        G_CALCs{VV} = G_calc;
        ZRADS{VV} = Zrad;
        Z1S{VV} = Z_1_calc;
        FREQS{VV} = fsweep;
        Vel{VV} = vel_calc;
        Fpeak{VV} = Fpeak_calc;
        xdef{VV} = x_calc;
        Zin{VV} = Z_calc;%Z_calc
        KAXX{VV} = kax;
        Cmeas{VV} = C_meas;
    end %if calc freq response
end %end voltage for loop

%% Extract output parameters

%Determine: center frequency fc = Max(Re(Zin))
%Then determine resonance frequency fr = Min(abs(Zin))
%Then determine antiresonance frequency fa = Max(abs(Zin))
%note: fr and fa are local minima. High low frequency impedance or low high frequency impedance can cause you to obtain incorrect fa and fr.
%solution: find fc, first. We know fr< fc < fa. So consider the range from 0-fc for fr and from fc-fmax for fa.

if calcFreqResponse
    for VV = 1:1:length(Vop)
        [~,FcIDX(VV)] = max(real(Zin{VV}));
        resFreq(VV) = FREQS{VV}(FcIDX(VV));
    end
    %Output fc as both "resFreq" and "fc" for backwards compatibility with old scripts.
    Sout.resFreq = resFreq;
    Sout.fc = resFreq;
    
    %get resonance freq fr
    for VV = 1:1:length(Vop)
        [~,FrIDX(VV)] = min(abs(Zin{VV}(1:FcIDX(VV))));
        ResonantFreq(VV) = FREQS{VV}(FrIDX(VV));
    end
    Sout.fr = ResonantFreq;

    %antiresonance freq: max(abs(Zin))
    for VV = 1:1:length(Vop)
        [~,FaIDX] = max(abs(Zin{VV}(FrIDX(VV):end))); %gives the index from fr:end
        FaIDX = FaIDX + FrIDX(VV)-1; %add FrIDX to FaIDX to get the correct index (subtrack 1 because 1 indexing)
        if FaIDX > length(FREQS{VV}) %stops weird edge case where fr is equal to fmax from giving error
           FaIDX = length(FREQS{VV});
        end
        AresFreq(VV) = FREQS{VV}(FaIDX);
    end
    Sout.fa = AresFreq;

    %Other method for obtaining resonant frequency: Maximum velocity
    for VV = 1:1:length(Vop)
        [~,maxIDX] = max(abs(Vel{VV}));
        ActualresFreq2(VV) = FREQS{VV}(maxIDX);
    end
    Sout.fr2 = ActualresFreq2;

    %Get average quantities and pressure:
    Ppeak = Fpeak;
    for VV = 1:length(Vop)
        %through variables for each domain:
        vRMS{VV} = Vel{VV} * alpha; %RMS average velocity
        vAvg{VV} = Vel{VV} * beta; %average velocity
        v0{VV} = Vel{VV}; %peak velocity;
        %across variables for each domain (Output force)
        F_RMS{VV} = Fpeak{VV} / alpha; %force variable for RMS domain
        F_avg{VV} = Fpeak{VV} / beta; %force variable for avg domain
        %F0{VV} = Fpeak{VV}; %force variable for peak domain

        %Pressure output:
        %from Koymen paper in references, particle velocity, force, and pressure in acoustic medium match RMS domain:
        %also Pout = F_RMS/A
        vout{VV} = vRMS{VV};
        Fout{VV} = F_RMS{VV};
        Pout{VV} = F_RMS{VV}/A;
        %average pressure and peak pressure may not have physical meaning but are sometimes useful to compare with simulations:
        P_avg{VV} = F_avg{VV}/A;
        Ppeak{VV} = Fpeak{VV} ./ (A);
        PowerOut{VV} = F_RMS{VV}.*vRMS{VV};
    end
    

    
    %Storing important Frequency outputs
    %length of fsweep varies by voltage - need to use cell arrays to return
    Sout.f = FREQS;
    Sout.Kax = KAXX;
    Sout.G = G_CALCs;
    Sout.G1 = GM_CALCs; %this is Just Z1
    Sout.ZRADS = ZRADS; %radiation impedance (legacy name)
    Sout.ZR0 = ZRADS; %radiation impedance
    Sout.Z1 = Z1S;
    Sout.Zin = Zin;
    Sout.Vel = Vel;
    Sout.xdef =xdef;
    Sout.w0ac =xdef;
    Sout.Fpeak = Fpeak;
    Sout.Cmeas = Cmeas; %"measured capacitance"
    %other parameters:with some duplicate outputs for convenience/backwards compatibility
    %velocity:
    Sout.vRMS = vRMS;
    Sout.vAvg = vAvg;
    Sout.v0 = v0;
    Sout.vOut = vout;
    %Force:
    Sout.Fout = Fout;
    Sout.FRMS = F_RMS;
    Sout.FAvg = F_avg;
    Sout.F0 = Fpeak;
    %Pressure:
    Sout.Pout = Pout;
    Sout.PRMS = Pout;
    Sout.PAvg = P_avg;
    Sout.Ppeak = Ppeak;
    Sout.P0 = Ppeak;
    Sout.P = Pout; %lumped output pressure
    Sout.PowerOut = PowerOut;
    Sout.fpred = fpred;
end

%Static deflection with increasing and decreasing voltage for hysteresis
    Sout.VINCLong = VINCLong;
    Sout.w0INCLong = w0INCLong;
    Sout.VDC = VINCLong;
    Sout.w0 = w0INCLong;
    Sout.C0 = CINCLong;
    Sout.CINCLong = CINCLong;
    Sout.CINC = CINC;
    Sout.VINC = VINC;
    Sout.w0INC = w0INC;
if calcHyst
    Sout.VDEC = VDEC;
    Sout.VDECLong = VDECLong;
    Sout.Vsb =Vsb;
    Sout.w0DEC = w0DEC;
    Sout.w0DECLong = w0DECLong;
    Sout.CDEC = CDEC;
    Sout.CDECLong = CDECLong;
end

Sout.Water = Water; %really this is pointless to return but it is convenient
if DeviceType > 0
    if strcmp(Sout.PostContactMode,'touchdown') %return the real collapse voltage and the touchdown voltage
        Sout.Vt = Vt;
        Sout.Vc = VcReal; %note: this collapse voltage is an over-estimate.
    else
        Sout.Vc = Vc;
    end
else
    Sout.Vc = Vc;
end

Sout.Params = CMUTParams;
Sout.wPI = wPI;

if CMUTParams.SimpleOutputs
    %remove legacy outputs that may confuse new users, add some new simple outputs.    
    if calcFreqResponse
        %remove extra frequency response outputs:
        Sout = rmfield(Sout,'resFreq'); %old name for fc
        Sout = rmfield(Sout,'xdef'); %old name for w0ac
        Sout = rmfield(Sout,'Vel'); %old name for v0
        Sout = rmfield(Sout,'Fpeak'); %old name for F0
        Sout = rmfield(Sout,'Ppeak'); %old name for P0
        % Sout = rmfield(Sout,'fr2'); %Probably not needed.
        Sout = rmfield(Sout,'ZRADS'); %Probably not needed.
        Sout = rmfield(Sout,'fpred');
        Sout = rmfield(Sout,'Cmeas');
        %Impedance of the mechanical portion:
        Sout = rmfield(Sout,'G1');
        Sout = rmfield(Sout,'Z1');   
    end

    if calcHyst
        Sout.VDEC = VDECLong;
        Sout.w0DEC = w0DECLong;
        Sout.CDEC = CDECLong;
        %Sout = rmfield(Sout,{'VINCLong','w0INCLong','VINC','w0INC','VDECLong','w0DECLong'});
        Sout = rmfield(Sout,{'VINCLong','w0INCLong','VDECLong','w0DECLong','CINCLong','CDECLong'});
    else
        Sout = rmfield(Sout,{'VINCLong','w0INCLong','CINCLong'});
    end

    if (DeviceType > 0) && strcmp(Sout.PostContactMode,'touchdown')
        % In the case of touchdown, this code also calculates a collapse voltage assuming the membrane does not deflect further.
        % However, this value can be significantly over-estimated and may confuse inexperienced users.
        % Hence, this value is not output, since it is likely inaccurate anyways.
        Sout.Vc = rmfield(Sout,'Vc'); 
    end
end

%% Make Basic Parameter structure to output a less overwhelming list of parameters for new users:
BasicParams = CMUTParams;
BasicParams = rmfield(BasicParams,{'SuppressPrints','PrintInputParams','WarningMSGBox','SimpleOutputs'}); %Remove UX parameters
%Remove analysis options with default values:
if calcHyst; BasicParams = rmfield(BasicParams,{'calcHysteresis'});  end 
if calcFreqResponse; BasicParams = rmfield(BasicParams,{'calcFreqResponse'});  end 
if ~PistonZ; BasicParams = rmfield(BasicParams,{'UsePistonRadiationImpedance'});  end 
if CMUTParams.Clamped1DRadiationImpedance == 1; BasicParams = rmfield(BasicParams,{'Clamped1DRadiationImpedance'});  end     
if CMUTParams.BiLayerSpringConstant; BasicParams = rmfield(BasicParams,{'BiLayerSpringConstant'});  end  
%remove un-used bias voltage info:
if BiasByCollapsePercent; BasicParams = rmfield(BasicParams,{'ManualBiasVoltage'});else;BasicParams = rmfield(BasicParams,{'collapsePercent'});end 
%Remove un-used & default frequency sweep settings:
if ~calcFreqResponse; BasicParams = rmfield(BasicParams,{'UseDefaultFrequencySweep','CustomFrequencySweepRange','fmin','fmax','fstep','fstepFine','EnableFineSweepInWater','fineFsweepRangeScalingFactor'});
else
    if CMUTParams.UseDefaultFrequencySweep; BasicParams = rmfield(BasicParams,{'UseDefaultFrequencySweep','CustomFrequencySweepRange'});  end
    if Water && ~CMUTParams.EnableFineSweepInWater; BasicParams = rmfield(BasicParams,{'fstepFine','EnableFineSweepInWater','fineFsweepRangeScalingFactor'});  end
    if ~Water; BasicParams = rmfield(BasicParams,{'EnableFineSweepInWater'});  end
end
  
%Remove un-used circuit parameters:
if ~calcFreqResponse; BasicParams = rmfield(BasicParams,{'VAC','EnableSeriesResistor','Rs','EnableSeriesCapacitor','Cs','EnableParasiticCapacitance','Cp'});
else
    if ~SeriesResistor && calcFreqResponse; BasicParams = rmfield(BasicParams,{'EnableSeriesResistor','Rs'});  end
    if ~CMUTParams.EnableSeriesCapacitor; BasicParams = rmfield(BasicParams,{'EnableSeriesCapacitor','Cs'});  end
    if ~CMUTParams.EnableParasiticCapacitance; BasicParams = rmfield(BasicParams,{'EnableParasiticCapacitance','Cp'});  end
end
%Remove un-used fluid medium parameters:
if ~calcFreqResponse; BasicParams = rmfield(BasicParams,{'Water','UseCustomFluidMedium','FluidDensity','FluidSpeedOfSound'});
elseif ~CMUTParams.UseCustomFluidMedium; BasicParams = rmfield(BasicParams,{'FluidDensity','FluidSpeedOfSound'});  end 
%Keep all core device parameters. Remove empty metal layers
if CMUTParams.h_metal==0; BasicParams = rmfield(BasicParams,{'rho_metal','E_metal'});  end 
if CMUTParams.h_metal2==0; BasicParams = rmfield(BasicParams,{'h_metal2','rho_metal2'});  end 
%Remove unused EP and IIP parameters: 
%CD:
if DeviceType == 0; BasicParams = rmfield(BasicParams,{'circularPostAreaCorrection','ap','d1','at','n_1postRows','n_2postRows','xp','SetPostsAtCollapseHeight'});  end 
if DeviceType == 1; BasicParams = rmfield(BasicParams,{'at'});  end %EP:
if (DeviceType > 0) && (CMUTParams.n_2postRows == 0) ; BasicParams = rmfield(BasicParams,{'xp'});  end %no 2 post rows = xp doesn't matter
%Remove un-used electrode trench modeling parameters:
if ~CMUTParams.ModelIsolatingTrench; BasicParams = rmfield(BasicParams,{'ModelIsolatingTrench','IsolatingTrenchWidth'});  end

Sout.BasicParams = BasicParams;

%% Output Sout

CMUTNums = Sout;%return all these outputs

end %end function