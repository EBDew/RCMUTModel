%Reproduce Fig 14 in our Paper:
%Add the path of the modeling code: (requires running the whole script)
filePath = matlab.desktop.editor.getActiveFilename; %get active file location
TutorialsFolder = fileparts(filePath); ProjectFolder = fullfile(TutorialsFolder,'..');
CodePath = fullfile(ProjectFolder,"Code/");addpath(genpath(CodePath))
FigureDataPath = fullfile(ProjectFolder,"Data/Figures/");addpath(genpath(FigureDataPath))

%% Regular Sized CMUT Model Parameters:
%define a struct with our device dimensions & analysis parameters
clear CMUT; CMUT = struct;

%device dimensions: 
CMUT.a = 53.5e-6; %membrane half-width [m]
CMUT.h = 4.9e-6; %membrane thickness [m]
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
CMUT.BiLayerSpringConstant = 1; % 1 (default) - accounts for electrode stiffness when calculating metal layer. 0 - neglects electrode layer stiffness

%Frequency Sweep Parameters:
%Using default frequency sweep:
CMUT.fmin = 100e3; %Minimum Frequency [Hz]
CMUT.fmax = 5e6; %Maximum Frequency [Hz]
CMUT.fstep = 5e3; %Frequency step size [Hz]
CMUT.fstepFine = 0.25e3; %Fine frequency step size [Hz] (air only).

%oil Params:
CMUT.Water = 1;
CMUT.UseCustomFluidMedium = 1;
CMUT.FluidDensity = 920.29; %COMSOL properties
CMUT.FluidSpeedOfSound = 1484.1; %This this one is from a paper @ 25 deg
%CMUT.SimpleOutputs = 0;

%Radiation impedance (Default)
CMUT.Clamped1DRadiationImpedance = 1; %1 - Use Radiation impedance corresponding to 1D velocity profile (recommended, default)
CMUT.VAC = 5; %Amplitude of driving sinusoidal signal [V]
%Define Operating points for the Frequency Sweep
CMUT.ManualBiasVoltage = [20 50 100 150 175];
CMUT.BiasByCollapsePercent = 0;
plotIdx = 2; %plot @ bias voltage #2, ie 50V

%% Large CMUT Model Parameters:
%Copy the parameters from CMUT, then change dimensions
CMUT2 = CMUT;
CMUT2.a = 250e-6;
CMUT2.h = 60e-6;
CMUT2.d0 = 275e-9;

%% Plot Settings:
plotPAbs = 1;
plotStatic = 1;
plotVRMS = 1;
plotPower = 1;

PlotANSYS_ES = 1; %Add FEM to static graph
PlotCOMSOL = 1; %Add FEM to dynamic graphs

plotFontSize = 16;
simpleYlabel = 0; %If enabled labels are all text which makes the figures the same size,

Plot3DDeflection = 1;

%% Simulate Regular Sized CMUT:
test = edrc(CMUT);
P1 = abs(test.P{plotIdx});
[P1Max, p1MaxIdx] = max(P1);
Fc1 = test.f{plotIdx}(p1MaxIdx);

%% Simulate Large Sized CMUT:
test2 = edrc(CMUT2);
P2 = abs(test2.P{plotIdx});
[P2Max, p2MaxIdx] = max(P2);
Fc2 = test2.f{plotIdx}(p2MaxIdx);

lambda_2 = CMUT.FluidSpeedOfSound / Fc2 ;

Frac_2 = 2*CMUT2.a / lambda_2;

fprintf('\nSize of Large device: %.1f (um)\n',2*CMUT2.a*1e6)
fprintf('Center Frequency: Theory, %.2f (MHz).\n',Fc2/1e6);
fprintf('Size WRT Lambda: %.2f.\n', Frac_2);
fprintf('Size WRT Lambda: lambda/%i.\n',round(1/Frac_2));
fprintf('Collapse voltage: %.1f V.\n',test2.Vc);
fprintf('Force Increase Ratio: %.2f .\n', max(abs(test2.Fout{plotIdx}))/max(abs(test.Fout{plotIdx}))   );
fprintf('Far-Field Pressure Increase Ratio: %.2f .\n', (P2Max * CMUT2.a) / ((P1Max * CMUT.a)));
fprintf('Acoustic Power Increase Ratio: %.2f .\n', max(abs(test2.Fout{plotIdx}.*test2.vRMS{plotIdx}))  / max(abs(test.Fout{plotIdx}.*test.vRMS{plotIdx}))  );
fprintf('Surface Pressure Ratio: %.2f %%.\n',P2Max/P1Max * 100);

%% Load FEM Data:
if PlotANSYS_ES
    load("Fig14ANSYS.mat")
end
if PlotCOMSOL
    load("Fig14COMSOL.mat")
end

%% Plotting:
ScreenSize = get(0,'ScreenSize');
FigWidth = 500;
FigHeight = 400;
XStartPos = 800;
YStartPos = 595;
NoFigs = 0;
figoffset = 222;

CustomColours(1,:) = 0.15* [1 1 1];
CustomColours(2,:) = 1/255 * [72 118 255]; %blue
CustomColours(3,:) = 1/255 * [220 20 60]; %red
CustomColours(4,:) = 1/255 * [255 215 0]*0.9; %Gold
CustomColours(5,:) = 1/255 * [139 0 139]; %purple
CustomColours(6,:) = 1/255 * [0 205 102]; %green

if plotPAbs
    LargeMemFig = figure(figoffset+NoFigs);
    clf;
    try
        LargeMemFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        LargeMemFig.Position = [(XStartPos + NoFigs*FigWidth) YStartPos FigWidth FigHeight];
    end   
    NoFigs = NoFigs + 1;
    hold on;box on;
    
    plot(test.f{plotIdx}/1e6,P1/1e3,'-','Color',CustomColours(1,:),'linewidth',2)
    plot(test2.f{plotIdx}/1e6,P2/1e3,'-','Color',CustomColours(2,:),'linewidth',2)
    if simpleYlabel
        ylabel('Pressure Magnitude (kPa)')
    else
        ylabel('Pressure |p_o_u_t| (kPa)')
    end    

    xlabel('Frequency (MHz)')

    legend('λ/8 device','λ/2 device')
    legend('boxoff');
    if PlotCOMSOL
        plot(f3_tot,abs(P3RMS)/1e3,':','Color',CustomColours(3,:),'linewidth',2)
        legend('λ/8 device (model)','λ/2 device (model)','λ/2 device (FEM)')
    end

    ax1 = gca;
    ax1.FontSize = plotFontSize;
    ax1.LineWidth = 1.75;
    ax1.FontName = 'Arial';
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    title('Surface Pressure');
end

if plotStatic
    %Collapse Voltage Plot:
    VcFig = figure(figoffset+NoFigs);
    clf;
    try
        VcFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        VcFig.Position = [(XStartPos + NoFigs*FigWidth) YStartPos FigWidth FigHeight];
    end
    NoFigs = NoFigs + 1;
    plot(test.VDC,test.w0*1e9,'-','Color',CustomColours(1,:),'linewidth',2)
    hold on;
    plot(test2.VDC,test2.w0*1e9,'-','Color',CustomColours(2,:),'linewidth',2)
    legend('λ/8 device','λ/2 device')
    if PlotANSYS_ES
        plot(VoltageRange,-MaxDisp,':sq','Color',CustomColours(3,:),'linewidth',2)        
        legend('λ/8 device (model)','λ/2 device (model)','λ/2 device (FEM)')
    end
    legend('boxoff');
    legend('Location','northwest')
    xlabel('Bias Voltage (V)')

    if simpleYlabel
        ylabel('Peak Deflection w0 (nm)')
    else
        ylabel('Peak Deflection w_0 (nm)')
    end
    
    ax2 = gca;
    ax2.FontSize = plotFontSize;
    ax2.LineWidth = 1.75;
    ax2.FontName = 'Arial';
    ax2.XColor = 'k';
    ax2.YColor = 'k';
    title('Electrostatic Behaviour');
    xlim([0 200])
end

if plotVRMS
    %Plot RMS Velocity
    VelFig = figure(figoffset+NoFigs);
    clf;
    try
        VelFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        VelFig.Position = [(XStartPos + (NoFigs-2)*FigWidth) (YStartPos - FigHeight - 85) FigWidth FigHeight];
    end    
    NoFigs = NoFigs + 1;
    plot(test.f{plotIdx}/1e6,abs(test.vRMS{plotIdx}*1e3),'-','Color',CustomColours(1,:),'linewidth',2)
    hold on;
    plot(test2.f{plotIdx}/1e6,abs(test2.vRMS{plotIdx}*1e3),'-','Color',CustomColours(2,:),'linewidth',2)
    xlabel('Frequency (MHz)')

    if simpleYlabel
        ylabel('Velocity Magnitude (mm/s)')
    else
        ylabel('Velocity |v_r_m_s| (mm/s)')
    end
    
    legend('λ/8 device','λ/2 device')
    if PlotCOMSOL
        plot(f3_tot,abs(v3RMS)*1e3,':','Color',CustomColours(3,:),'linewidth',2)
        legend('λ/8 device (model)','λ/2 device (model)','λ/2 device (FEM)')
    end

    legend('boxoff');
    ax4 = gca;
    ax4.FontSize = plotFontSize;
    ax4.LineWidth = 1.75;
    ax4.FontName = 'Arial';
    ax4.XColor = 'k';
    ax4.YColor = 'k';
    title('RMS Membrane Velocity')
end

if plotPower
    %Plot Power:
    PFig = figure(figoffset+NoFigs);
    clf;
    try
        PFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        PFig.Position = [(XStartPos + (NoFigs-2)*FigWidth) (YStartPos - FigHeight - 85) FigWidth FigHeight];
    end
    NoFigs = NoFigs + 1;
    plot(test.f{plotIdx}/1e6,abs(test.vRMS{plotIdx}.*test.FRMS{plotIdx})*1e3,'-','Color',CustomColours(1,:),'linewidth',2)
    hold on;
    plot(test2.f{plotIdx}/1e6,abs(test2.vRMS{plotIdx}.*test2.FRMS{plotIdx})*1e3,'-','Color',CustomColours(2,:),'linewidth',2)
    xlabel('Frequency (MHz)')

    if simpleYlabel
        ylabel('Power Magnitude (mW)')
    else
        ylabel('Power |F_r_m_s \times v_r_m_s| (mW)')
    end  

    legend('λ/8 device','λ/2 device')
    if PlotCOMSOL
        plot(f3_tot,abs(v3RMS .* F3RMS)*1e3,':','Color',CustomColours(3,:),'linewidth',2)
        legend('λ/8 device (model)','λ/2 device (model)','λ/2 device (FEM)')
    end
    legend('boxoff');
    ax3 = gca;
    ax3.FontSize = plotFontSize;
    ax3.LineWidth = 1.75;
    ax3.FontName = 'Arial';
    ax3.XColor = 'k';
    ax3.YColor = 'k';
    title('Acoustic Power Output')
end

if Plot3DDeflection && PlotANSYS_ES
    %Reconstruct 3D profile from X and Y profile
    ii = 6;
    XZNorm{ii} = -XZValues{ii} ./max(-XZValues{ii});
    YZNorm{ii} = -YZValues{ii} ./max(-YZValues{ii});
    XYZNorm{ii} = XZNorm{ii} .* YZNorm{ii}'; %120x200 grid for x and Y
    
    %plot settings:
    shrinkfactor = 1; %This is useful for plotting larger aspect ratios    
    FEMXYFig = figure(figoffset+NoFigs);
    clf;
    try
        FEMXYFig.Position = GenerateFigPositions(0,560,420, 100,300,ScreenSize,85);
    catch
        FEMXYFig.Position = [100 200 560 420];
    end
    NoFigs = NoFigs + 1;
    hold on
    
    imagesc(XValues{ii},YValues{ii}/shrinkfactor,XYZNorm{ii}'*1e2);
    xlabel('x (µm)');
    ylabel('y (µm)');
    axis image
    yticks([-1500 -1000 -500 0 500 1000 1500]/shrinkfactor);
    %define axes:
    ymin = -max([-min(YValues{ii}) max(YValues{ii})]);
    ymax = -ymin;
    xmin = -max([-min(XValues{ii}) max(XValues{ii})]);
    xmax = -xmin;
    yTix = [ymin:(ymax-ymin)/6:max(YValues{ii})]/shrinkfactor;
    yticks(yTix)    
    yticklabels(num2str(yTix'*shrinkfactor))
    xticks([xmin 0 xmax])

    axis image
    box on

    title('FEM')
    colorbar;
    ax1 = gca;
    ax1.FontSize = 14;
    ax1.FontName = 'Arial';
    ax1.XColor = 'k';
    ax1.YColor = 'k';

    colormap parula;
    set(ax1,"TickDir","out")
    set(ax1,'XTickLabelRotation',0)   

end