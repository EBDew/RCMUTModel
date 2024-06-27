%Template script for simulating an IIP rectangular CMUT (Electrode Posts):

%Add the path of the modeling code: (requires running the whole script)
filePath = matlab.desktop.editor.getActiveFilename; %get active file location
TutorialsFolder = fileparts(filePath); ProjectFolder = fullfile(TutorialsFolder,'..');
CodePath = fullfile(ProjectFolder,"Code/");addpath(genpath(CodePath))

%define a struct with our device dimensions & analysis parameters
clear CMUT; CMUT = struct;

%% Device Parameters

%device dimensions: 
CMUT.a = 60e-6; %membrane half-width [m]
CMUT.h = 5e-6; %membrane thickness [m]
CMUT.d0 = 480e-9; %cavity vacuum gap [m]
CMUT.b = 1.5e-3; %membrane half-length [m]
CMUT.di = 360e-9; %insulating layer thickness [m]

%material properties:
CMUT.epsr = 3.9; %relative permittivity of dielectric layer
CMUT.E = 148e9; %Young's Modulus of membrane [Pa]
CMUT.rho = 2329; %Density of membrane [kg/m^3]

%Metal electrode layer:
CMUT.h_metal = 250e-9; %Electrode thickness [m]
CMUT.E_metal = 76e9; %Young's Modulus of metal electrode [Pa]
CMUT.rho_metal = 19300; %Density of metal electrode [kg/m^3]

%IIP Specific Parameters:
CMUT.DeviceType = 'IIP'; %Indicate CMUT architecture. Default "CD" = no posts. EP = Electrode posts. IIP = Isolated Isolation Posts.
CMUT.ap  = 5e-6; %post radius [m]. Only used with DeviceType = "EP" or "IIP" 
CMUT.d1 = 180e-9; %distance to posts from undeflected membrane [m]. (EP or IIP only)
CMUT.at  = 2e-6; %trench width around each IIPs [m]. IIPs only.
CMUT.n_1postRows = 58; %number of rows with 1 central post (EP or IIP only)
CMUT.n_2postRows = 59; %number of rows with 2 posts. Each post is a distance of ± xp from the center. (EP or IIP only)
CMUT.xp = 25e-6; %location of center of the posts in the 2post rows, with respect to the center of the membrane, the post centers are at ±xp. (EP or IIP only)

CMUT.SetPostsAtCollapseHeight = 0; %Runs a CD static simulation to determine the collapse point, then sets d1 equal to that pullin deflection. Overrules d1 value. (EP or IIP only)



%% Dynamic Analysis Parameters

%Acoustic medium
CMUT.Water = 1; %Toggle between air or water. 1 = Water, 0 = air. Default - 0 (air).
CMUT.UseCustomFluidMedium = 1; %1 = enabled, 0 = disabled. Default - Disabled.
% FluidDensity and FluidSpeedOfSound do not do anything unless UseCustomFluidMedium = 1.
CMUT.FluidDensity = 920.29; %Custom fluid medium density [kg/m^3]
CMUT.FluidSpeedOfSound = 1484.1; %Custom fluid medium speed of sound [m/s]

%Frequency Sweep Parameters:
%Using default frequency sweep:
CMUT.fmin = 100e3; %Minimum Frequency [Hz]
CMUT.fmax = 5e6; %Maximum Frequency [Hz]
CMUT.fstep = 10e3; %Frequency step size [Hz]
CMUT.fstepFine = 2e3; %Fine frequency step size [Hz] (air only).
%This will define a frequency array from fmin to fmax in steps of fstep. 
%In the case of air (CMUT.Water = 0), points near the predicted resonance frequency will
%have a finer frequency step size determined by CMUT.fstepFine.

%Enable & define custom frequency sweep: (uncomment below lines)
% CMUT.UseDefaultFrequencySweep = 0; % 1 = enable default sweep, 0 = use custom sweep specified by CustomFrequencySweepRange
% CMUT.CustomFrequencySweepRange = 100e3:5e3:5e6; %Custom Frequency Sweep. Enabled with UseDefaultFrequencySweep = 0

%Define Operating points for the Frequency Sweep: 2 Options:
%option 1: Specify fractions of the collapse voltage:
CMUT.BiasByCollapsePercent = 1; % 1 - Uses collapsePercent (default), 0 - uses ManualBiasVoltage
CMUT.collapsePercent = [0.2 0.5 0.8 0.95]; %Specify bias voltage as fraction of collapse voltage. Must be greater than 0 and less than 1.
%option 2: (Uncomment below). Specify bias voltages:
% CMUT.BiasByCollapsePercent = 0; %uses ManualBiasVoltage
% CMUT.ManualBiasVoltage = 20:20:120;%Operating voltage for dynamic model in [V], must be less than collapse voltage

CMUT.VAC = 1; %Amplitude of driving sinusoidal signal [V]
%% Advanced Parameters:
% It is usually best to leave these at their default values.

% Simple outputs:
CMUT.SimpleOutputs = 1; %intended for new users. Hides some extra outputs that may be confusing.
% Most of these outputs are legacy outputs from old versions of the code

%Frequency sweep extra features:
CMUT.EnableFineSweepInWater = 0; %set this to 1 to enable the fine sweep mechanics in water. 
%(EnableFineSweepInWater is not recommended because it is difficult to predict the operating frequency before running the simulation.)
CMUT.fineFsweepRangeScalingFactor = 1; %Factor which the width of the fine sweep is multiplied by. Increase = more range swept with fstepFine.

% Radiation impedance:
CMUT.Clamped1DRadiationImpedance = 1; %1 - Use Radiation impedance corresponding to 1D velocity profile (recommended, default)
% [0 or 2] - Use radiation impedance corresponding to 2D square-like velocity profile (maybe appropriate at low aspect ratios)
CMUT.UsePistonRadiationImpedance = 0; % 1 - Enable, uses rectangular piston radiation impedance instead. Not recommended.

% Account for electrode stiffness
CMUT.BiLayerSpringConstant = 1; % 1 (default) - accounts for electrode stiffness when calculating metal layer. 0 - neglects electrode layer stiffness

%Model an isolating trench around bottom electrode as described in https://ieeexplore.ieee.org/document/10026669
CMUT.ModelIsolatingTrench = 0; %0 (default) disabled. 1 - enable. Makes a trench of width 'IsolatingTrenchWidth' around the bottom electrode 
CMUT.IsolatingTrenchWidth = 5e-6; %isolating width trench in [m]. Not enabled unless ModelIsolatingTrench = 1.

%Simulate parasitic circuit elements:
CMUT.EnableSeriesResistor = 0; %Model a resistor in series with the voltage source in the electrical domain (1 - enable, 0 disable)
CMUT.Rs = 38; %Series resistance value [ohms]. Unused unless EnableSeriesResistor = 1.
CMUT.EnableSeriesCapacitor = 0; %Model a capacitor in series with the voltage source in the electrical domain (1 - enable, 0 disable)
CMUT.Cs = 10e-9; %Series capacitance [Farads]. Unused unless EnableSeriesCapacitor = 1.
CMUT.EnableParasiticCapacitance = 0; %Model a parasitic capacitor in parallel with the voltage source in the electrical domain (1 - enable, 0 disable)
CMUT.Cp = 20e-12; %Parallel (parasitic) capacitance [Farads]. Unused unless EnableParasiticCapacitance = 1.

% EP/ IIP Only:
CMUT.circularPostAreaCorrection = 1; %1 - enable (default), 0 - disable. Scales the area of the modeled posts to be equivalent to the area of circular posts
% Enable: models posts as rectangles the same area as the circles., x dimensions: 2*ap, y-dimension: 2*ap*pi/4 (recommended)
% Disable :models posts as squares with the side length 2*ap

%% Simulate the device with the chosen parameters:
results = edrc(CMUT); %simulates the device with the chosen parameters, returns a structure with the results.
%% Print useful Metrics:
%EP/IIP only:
fprintf('\nPosts Area Percentage = %.1f %%',results.PostAreaFraction*100)
%Electrostatics:
fprintf(strcat('\nContact regime between membrane and posts:',results.PostContactMode,'.'))
if strcmp(results.PostContactMode,'touchdown')
    fprintf('\nTouchdown Voltage = %.1f V',results.Vt)
else %collapse    
    fprintf('\nCollapse Voltage = %.1f V',results.Vc)
end

fprintf('\nDeflection @ V = 0: %.1f nm',results.w0(1)*1e9)
fprintf('\nCapacitance @ V = 0: %.2f pF\n',results.C0(1)*1e12)

for ii = 1:length(results.Vop)
    fprintf('\nFrequency response @ V = %.1f V:',results.Vop(ii))
    fprintf('\nResonant Frequency fr = %.2f MHz',results.fr(ii)/1e6)
    fprintf('\nAntiresonance Frequency fa = %.2f MHz',results.fa(ii)/1e6)
end

%% Plot Outputs:
ScreenSizeFraction = 0.8; %fraction of the screen occupied by the results figure
screenSize = get(0,'ScreenSize');

ResultPlot = figure;
ResultPlot.Position = [screenSize(1)+20 screenSize(1)+50  screenSize(3:4)*ScreenSizeFraction];
T = tiledlayout('flow');

%Plot electrostatic analysis
nexttile
plot(results.VDC,results.w0*1e9,'k','LineWidth',1.5)
box on;
xlabel('Bias Voltage (V)')
ylabel('Membrane Deflection w_0 (nm)')
title('Electrostatic Analysis')

% get a string representing the acoustic medium:
if results.Params.UseCustomFluidMedium
    mediumString = 'Custom';
elseif results.Params.Water
    mediumString = 'Water';
else
    mediumString = 'Air';
end


%Plot Impedance magnitude
nexttile
xlabel('Frequency (MHz)')
ylabel('Input Impedance | Z_i_n | (ohms)')
hold on; box on;
legendstr = cell(length(results.Vop),1);
for ii = 1:length(results.Vop)
    plot(results.f{ii}/1e6,abs(results.Zin{ii}));
    legendstr{ii} = sprintf('Bias Voltage = %.1f V',results.Vop(ii));
end
set(gca, 'YScale', 'log')
legend(legendstr);legend('location','best');legend('boxoff')
% title('Input Impedance Magnitude')
title(strcat('Input Impedance Magnitude (',mediumString,' Acoustic Medium)'))



%Plot velocity magnitude
nexttile
xlabel('Frequency (MHz)')
ylabel('Membrane Velocity | v_0 | (mm/s)')
hold on; box on;
legendstr = cell(length(results.Vop),1);
for ii = 1:length(results.Vop)
    plot(results.f{ii}/1e6,abs(results.v0{ii})*1e3);
    legendstr{ii} = sprintf('Bias Voltage = %.1f V',results.Vop(ii));
end
legend(legendstr);legend('location','best');legend('boxoff')
% title('Membrane Velocity Magnitude')
title(strcat('Membrane Velocity Magnitude (',mediumString,' Acoustic Medium)'))

%Plot pressure magnitude
nexttile
xlabel('Frequency (MHz)')
ylabel('Pressure Output | P_o_u_t | (kPa)')
hold on; box on;
legendstr = cell(length(results.Vop),1);
for ii = 1:length(results.Vop)
    plot(results.f{ii}/1e6,abs(results.Pout{ii})/1e3);
    legendstr{ii} = sprintf('Bias Voltage = %.1f V',results.Vop(ii));
end
legend(legendstr);legend('location','best');legend('boxoff')
% title('RMS Output Pressure Magnitude')
title(strcat('RMS Output Pressure Magnitude (',mediumString,' Acoustic Medium)'))




