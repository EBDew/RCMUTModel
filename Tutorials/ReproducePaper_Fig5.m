%Reproduce Fig 5 in our Paper:

%Add the path of the modeling code: (requires running the whole script)
filePath = matlab.desktop.editor.getActiveFilename; %get active file location
TutorialsFolder = fileparts(filePath); ProjectFolder = fullfile(TutorialsFolder,'..');
CodePath = fullfile(ProjectFolder,"Code/");addpath(genpath(CodePath))
FigureDataPath = fullfile(ProjectFolder,"Data/Figures/");addpath(genpath(FigureDataPath))

%define a struct with our device dimensions & analysis parameters
clear CMUT; CMUT = struct;
%% model parameters:
%device dimensions: 
CMUT.a = 63.5e-6; %membrane half-width [m]
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
CMUT.Water = 0;

%Radiation impedance (Default)
CMUT.Clamped1DRadiationImpedance = 1; %1 - Use Radiation impedance corresponding to 1D velocity profile (recommended, default)
CMUT.VAC = 5; %Amplitude of driving sinusoidal signal [V]
%Define Operating points for the Frequency Sweep
CMUT.ManualBiasVoltage = [20 90];
CMUT.BiasByCollapsePercent = 0;

test = edrc(CMUT);

load Fig5_FEMData.mat %load FEM data - may require a path.

%% Plotting Parameters:
plotIndices = [1 2];
% Impedance:
plotZMag = 1;
Zplot_LogScale = 1;
Plot3DFigs = 1;


% Insets:
plotZMagInset2ndIdea = 1;
insetXlim = [2.6 3.0]; %Horizontal long part
shrinkfactor = 4.5;%5; %matches above
insetLW = 1.75; %matches other figure;

%% Plotting:
ScreenSize = get(0,'ScreenSize');
FigWidth = 500;
FigHeight = 400;
XStartPos = 200;
YStartPos = 595;
NoFigs = 0;
figoffset = 800;

CustomColors = 1;
%CustomColours(1,:) = 1/255 * [0 0 0];
CustomColours(1,:) = 0.15* [1 1 1];
CustomColours(2,:) = 1/255 * [72 118 255]; %blue
CustomColours(3,:) = 1/255 * [220 20 60]; %red
CustomColours(4,:) = 1/255 * [255 215 0]*0.9; %Gold
CustomColours(5,:) = 1/255 * [139 0 139]; %purple
CustomColours(6,:) = 1/255 * [0 205 102]; %green

clear LegStrings;clear Model3Colors;clear COMSOL3Colors;

%plot parameters:
C3LW = 1.25;
%test:
Marker2D = '-';
SwapFEMorder = 0;
    
COMSOLLW = [1.5; 1.5];


Model3Colors(1,:) = CustomColours(5,:);
Model3Colors(2,:) = CustomColours(2,:);

COMSOL3Colors(2,:) = 1/0.9*CustomColours(3,:);
COMSOL3Colors(1,:) = CustomColours(4,:);

COMSOLColors(2,:) = 0.5*CustomColours(1,:);
COMSOLColors(1,:) = 3* CustomColours(1,:);

%test7:
Marker3D = '--';
ModelLineW = 2.75;

if plotZMag
    clear LegStrings;
    ZmagFig = figure(figoffset+NoFigs);
    clf;
    try
        ZmagFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        ZmagFig.Position = [(XStartPos + NoFigs*FigWidth) YStartPos FigWidth FigHeight];
    end
    NoFigs = NoFigs + 1;
    hold on;box on;
    %for ii = 1:length(plotIndices) 
    for ii = length(plotIndices):-1:1
        if Plot3DFigs %include the 3D fig
            plot(test.f{plotIndices(ii)}/1e6,abs(test.Zin{plotIndices(ii)}),'-','Color',Model3Colors(ii,:),'linewidth',ModelLineW)
            plot(f_tot(:,plotIndices(ii)),ZMag(:,plotIndices(ii)),Marker2D,'Color',COMSOLColors(ii,:),'linewidth',COMSOLLW(ii));
            plot(f3_tot(:,plotIndices(ii)),Z3Mag(:,plotIndices(ii)),Marker3D,'Color',COMSOL3Colors(ii,:),'linewidth',C3LW)

            %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
            LegStrings{3*(ii-1) +1} = sprintf('COMSOL 3D %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{3*(ii-1) +2} = sprintf('COMSOL 2D %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{3*(ii-1) +3} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
        else
            %order swapped:
            plot(test.f{plotIndices(ii)}/1e6,abs(test.Zin{plotIndices(ii)}),'-','Color',CustomColours(1+ii,:),'linewidth',2)
            plot(f_tot(:,plotIndices(ii)),ZMag(:,plotIndices(ii)),'-','Color',COMSOLColors(ii,:),'linewidth',COMSOLLW(ii));
            %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
            LegStrings{2*(ii-1) +1} = sprintf('COMSOL %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{2*(ii-1) +2} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
        end
    end
    %ylabel('Impedance Magnitude (Ohms)')
    ylabel('Impedance Magnitude (Ω)')
    xlabel('Frequency (MHz)')
    
    
    %if reversed:
    LegStrings = LegStrings(length(LegStrings):-1:1);

    legend(LegStrings); 
    legend('boxoff');

    ax1 = gca;
    ax1.FontSize = 14;
    ax1.LineWidth = 1.75;
    %ax1.Legend.FontSize = 14;
    ax1.FontName = 'Arial';
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    title('Model vs FEM (Air)');

    if Zplot_LogScale
        set(gca, 'YScale', 'log')
    end
end


if plotZMag && plotZMagInset2ndIdea
   clear LegStrings;
    ZmagFigI = figure(figoffset+NoFigs);
    clf;
    try
        ZmagFigI.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        ZmagFigI.Position = [(XStartPos + NoFigs*FigWidth) YStartPos FigWidth FigHeight];
    end
    NoFigs = NoFigs + 1;
    hold on;box on;
    %for ii = 1:length(plotIndices) 
    for ii = length(plotIndices):-1:1
        if Plot3DFigs %include the 3D fig
            %order swapped:
            %this plot looks better with the 3D in front
            plot(test.f{plotIndices(ii)}/1e6,abs(test.Zin{plotIndices(ii)}),'-','Color',Model3Colors(ii,:),'linewidth',ModelLineW*2)    
            plot(f_tot(:,plotIndices(ii)),ZMag(:,plotIndices(ii)),'-','Color',COMSOLColors(ii,:),'linewidth',COMSOLLW(ii)*2); 
            plot(f3_tot(:,plotIndices(ii)),Z3Mag(:,plotIndices(ii)),Marker3D,'Color',COMSOL3Colors(ii,:),'linewidth',C3LW*2)   
            %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
            LegStrings{3*(ii-1) +1} = sprintf('COMSOL 3D %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{3*(ii-1) +2} = sprintf('COMSOL 2D %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{3*(ii-1) +3} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
        else
            %order swapped:
            plot(test.f{plotIndices(ii)}/1e6,abs(test.Zin{plotIndices(ii)}),'-','Color',CustomColours(1+ii,:),'linewidth',2)
            plot(f_tot(:,plotIndices(ii)),ZMag(:,plotIndices(ii)),'-','Color',COMSOLColors(ii,:),'linewidth',COMSOLLW(ii));
            %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
            LegStrings{2*(ii-1) +1} = sprintf('COMSOL %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{2*(ii-1) +2} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
        end
    end
    ax2 = gca;
    ax2.FontSize = 14;
    ax2.LineWidth = 1.75;
    %ax1.Legend.FontSize = 14;
    ax2.FontName = 'Arial';
    ax2.XColor = 'k';
    ax2.YColor = 'k';

    ax2.LineWidth = insetLW;


    if Zplot_LogScale
        set(gca, 'YScale', 'log')
    end
    
    %inset settings:
    xlim(insetXlim);

    %remove ticks:
    set(gca,'XTick',[])
    set(gca,'YTick',[])

    BoxX = ax2.XLim;
    BoxY = ax2.YLim;

    %draw a box on the original fig:
    figure(ZmagFig)
    ax1.Legend.AutoUpdate='off';
    boxLw = 1.25;
    plot(BoxX,BoxY(1)*[1 1],'k','linewidth',boxLw);
    plot(BoxX(2)*[1 1],BoxY,'k','linewidth',boxLw);
    plot(BoxX(1)*[1 1],BoxY,'k','linewidth',boxLw);
    plot(BoxX,BoxY(2)*[1 1],'k','linewidth',boxLw);

    ZFigX = ax1.XLim;
    ZFigY = ax1.YLim;    
    %make ax2 the same aspect ratio
    ax2.Position = [ax2.Position(1:2)  (diff(BoxX)/diff(ZFigX))/(diff(log10(BoxY))/diff(log10(ZFigY)))  * ax2.Position(4)   ax2.Position(4)];
    %then shrink ax2 by a factor of 'shrink factor' (defined above):
    ax2.Position = [ax2.Position(1:2)  ax2.Position(3)/shrinkfactor   ax2.Position(4)/shrinkfactor];
    moveFactor = 0.15;
    ax2.Position = [ax2.Position(1) + moveFactor  ax2.Position(1) + moveFactor  ax2.Position(3)   ax2.Position(4)];

end
