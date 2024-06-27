%Reproduce Fig 6 in our Paper:
%Add the path of the modeling code: (requires running the whole script)
filePath = matlab.desktop.editor.getActiveFilename; %get active file location
TutorialsFolder = fileparts(filePath); ProjectFolder = fullfile(TutorialsFolder,'..');
CodePath = fullfile(ProjectFolder,"Code/");addpath(genpath(CodePath))
FigureDataPath = fullfile(ProjectFolder,"Data/Figures/");addpath(genpath(FigureDataPath))

%define a struct with our device dimensions & analysis parameters
clear CMUT; CMUT = struct;
%% model parameters:
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
CMUT.ManualBiasVoltage = [20 80 170];
CMUT.BiasByCollapsePercent = 0;

test = edrc(CMUT);

load Fig6FEMData.mat %load FEM data - may require a path.

%% Plotting Parameters:
plotIndices = [1 3]; %for first figure
plotConductance = 1;
plotGvsPistonZ = 1;
pistonPlotIndex = 1;
Plot3DFigs = 1;
Zplot_LogScale = 1;

%% Plotting:
ScreenSize = get(0,'ScreenSize');
FigWidth = 500;
FigHeight = 400;
XStartPos = 200;
YStartPos = 595;
NoFigs = 0;
figoffset = 810;

CustomColors = 1;
%CustomColours(1,:) = 1/255 * [0 0 0];
CustomColours(1,:) = 0.15* [1 1 1];
CustomColours(2,:) = 1/255 * [72 118 255]; %blue
CustomColours(3,:) = 1/255 * [220 20 60]; %red
CustomColours(4,:) = 1/255 * [255 215 0]*0.9; %Gold
CustomColours(5,:) = 1/255 * [139 0 139]; %purple
CustomColours(6,:) = 1/255 * [0 205 102]; %green

clear LegStrings;clear Model3Colors;clear COMSOL3Colors; 
%initial formatting ideas
COMSOLColors(1,:) = CustomColours(1,:);
COMSOLColors(2,:) = 2* CustomColours(1,:);
COMSOLColors(3,:) = 4* CustomColours(1,:);
COMSOLColors(4,:) = 6* CustomColours(1,:);

COMSOLLW = [1.75; 1.5; 1.25; 1.0];

Model3Colors(2,:) = CustomColours(2,:);
Model3Colors(1,:) = 0.66* CustomColours(2,:);
Model3Colors(3,:) = 1.25*1/255 * [72 118 0] + [0 0 1];
Model3Colors(4,:) = 1.5*1/255 * [72 118 0] + [0 0 1];

COMSOL3Colors(2,:) = CustomColours(3,:);
COMSOL3Colors(1,:) = 0.5*CustomColours(3,:);
COMSOL3Colors(3,:) = 1.25* 1/255 * [0 20 60] + [1 0 0];
COMSOL3Colors(4,:) = 1.5* 1/255 * [0 20 60] + [1 0 0];

C3LW = 1.25;
%test:
Marker3D = '-.';
ModelLineW = 2.5;
COffset = 0;
Marker2D = '-';
clear minYlim; clear maxYlim;
% a bunch of formatting options I tested out
if length(plotIndices) == 2
    COMSOL3Colors(1,:) = CustomColours(4,:);    
    %test:
    Model3Colors(1,:) = CustomColours(5,:);
    Model3Colors(2,:) = CustomColours(2,:); %return to original color
    COMSOLColors(2,:) = 1.5* CustomColours(1,:);
    COMSOLColors(2,:) = CustomColours(1,:);
    COMSOLColors(1,:) = 3* CustomColours(1,:);
    COMSOL3Colors(2,:) = 1/0.9*CustomColours(3,:);
    COMSOLColors(2,:) = 0.5*CustomColours(1,:);
    
    ModelLineW = 2.75;

    Marker3D = '-';
    SwapFEMorder = 1;
    COMSOLLW = [2; 2.; 1.25; 1.0]; %dashes too spaced with 2.25
    Marker2D = '--';
    C3LW = 1.5;
    COMSOLColors(1,:) = 1.25 * CustomColours(1,:);



end
if length(plotIndices) == 1
    COMSOL3Colors(1,:) = CustomColours(4,:);    
    Model3Colors(1,:) = CustomColours(5,:);
    Model3Colors(2,:) = CustomColours(2,:); %return to original color
    COMSOLColors(2,:) = 1.5* CustomColours(1,:);
    COMSOLColors(2,:) = CustomColours(1,:);
    COMSOLColors(1,:) = 3* CustomColours(1,:);
    COMSOL3Colors(2,:) = 1/0.9*CustomColours(3,:);
    COMSOLColors(2,:) = 0.5*CustomColours(1,:);
    Marker3D = '-';
    ModelLineW = 2.75;
    COffset = 1;
    minYlim = 1e-4;
    SwapFEMorder = 1;
    Marker2D = '--';
    COMSOLLW = [2; 1.5; 1.25; 1.0];
    C3LW = 1.5;
end

if plotConductance
    clear LegStrings;
    CondFig = figure(figoffset+NoFigs);
    clf;
    try
        CondFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        CondFig.Position = [(XStartPos + NoFigs*FigWidth) YStartPos FigWidth FigHeight];
    end
    NoFigs = NoFigs + 1;
    hold on;box on;
    %for ii = 1:length(plotIndices) 
    for ii = length(plotIndices):-1:1
        if Plot3DFigs %include the 3D fig
            %order swapped:
            %this plot looks better with the 3D in front
            
            if SwapFEMorder %2d dashed on top
                %plot(test.f{plotIndices(ii)}/1e6,abs(test.xdef{plotIndices(ii)})*1e9,'-','Color',CustomColours(1+ii,:),'linewidth',1.75)
                plot(test.f{plotIndices(ii)}/1e6,real(test.G{plotIndices(ii)})*1e6,'-','Color',Model3Colors(ii+COffset,:),'linewidth',ModelLineW)
                plot(f3_tot(:,plotIndices(ii)),G3(:,plotIndices(ii))*1e6,Marker3D,'Color',COMSOL3Colors(ii+COffset,:),'linewidth',C3LW)
                plot(f_tot(:,plotIndices(ii)),G2(:,plotIndices(ii))*1e6,Marker2D,'Color',COMSOLColors(ii+COffset,:),'linewidth',COMSOLLW(ii));                
                %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
                LegStrings{3*(ii-1) +2} = sprintf('COMSOL 3D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{3*(ii-1) +1} = sprintf('COMSOL 2D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{3*(ii-1) +3} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
            else
                %plot(test.f{plotIndices(ii)}/1e6,abs(test.xdef{plotIndices(ii)})*1e9,'-','Color',CustomColours(1+ii,:),'linewidth',1.75)
                plot(test.f{plotIndices(ii)}/1e6,real(test.G{plotIndices(ii)})*1e6,'-','Color',Model3Colors(ii+COffset,:),'linewidth',ModelLineW)
                plot(f_tot(:,plotIndices(ii)),G2(:,plotIndices(ii))*1e6,Marker2D,'Color',COMSOLColors(ii+COffset,:),'linewidth',COMSOLLW(ii));
                plot(f3_tot(:,plotIndices(ii)),G3(:,plotIndices(ii))*1e6,Marker3D,'Color',COMSOL3Colors(ii+COffset,:),'linewidth',C3LW)
                %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
                LegStrings{3*(ii-1) +1} = sprintf('COMSOL 3D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{3*(ii-1) +2} = sprintf('COMSOL 2D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{3*(ii-1) +3} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
            end





        else
            %order swapped:
            plot(test.f{plotIndices(ii)}/1e6,real(test.G{plotIndices(ii)})*1e6,'-','Color',CustomColours(1+ii,:),'linewidth',2)
            plot(f_tot(:,plotIndices(ii)),G2(:,plotIndices(ii))*1e6,'-','Color',COMSOLColors(ii,:),'linewidth',COMSOLLW(ii));
            %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
            LegStrings{2*(ii-1) +1} = sprintf('COMSOL %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{2*(ii-1) +2} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));

        end
    end
    %ylabel('Conductance (Ohms^-^1)')
    ylabel('Conductance (µΩ^-^1)')
    
    xlabel('Frequency (MHz)')
    
    %if reversed:
    LegStrings = LegStrings(length(LegStrings):-1:1);

    legend(LegStrings); 
    legend('boxoff');
    legend('location','southeast');

    ax1 = gca;
    ax1.FontSize = 14;
    ax1.LineWidth = 1.75;
    %ax1.Legend.FontSize = 14;
    ax1.FontName = 'Arial';
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    title('Model vs FEM (Oil)');
    
    if Zplot_LogScale
        set(gca, 'YScale', 'log')
    end 
    
    if exist("minYlim")
        if exist("maxYlim")
            ylim([minYlim maxYlim])
        else
            ylim([minYlim inf])
        end
    elseif exist("maxYlim")
        ylim([-inf maxYlim])
    end
    

end

if plotConductance && plotGvsPistonZ
    pistonCMUT = CMUT;
    pistonCMUT.fmin=145e3; %avoids a bug where the program extrapolates the piston data past the lowest pre-calculated frequency
    pistonCMUT.PistonZ = 1;

    testP = edrc(pistonCMUT);
    plotIndices = pistonPlotIndex;
    COffset = 1;

    clear LegStrings;
    GFig = figure(figoffset+NoFigs);
    clf;
    try
        GFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
    catch
        GFig.Position = [(XStartPos + NoFigs*FigWidth) YStartPos FigWidth FigHeight];
    end  
    NoFigs = NoFigs + 1;
    hold on;box on;
    %for ii = 1:length(plotIndices) 
    for ii = length(plotIndices):-1:1
        if Plot3DFigs %include the 3D fig
            %order swapped:
            %this plot looks better with the 3D in front
            

            if SwapFEMorder %2d dashed on top
                %plot(test.f{plotIndices(ii)}/1e6,abs(test.xdef{plotIndices(ii)})*1e9,'-','Color',CustomColours(1+ii,:),'linewidth',1.75)
                plot(test.f{plotIndices(ii)}/1e6,real(test.G{plotIndices(ii)})*1e6,'-','Color',Model3Colors(ii+COffset,:),'linewidth',ModelLineW)
                plot(testP.f{plotIndices(ii)}/1e6,real(testP.G{plotIndices(ii)})*1e6,'-','Color',CustomColours(4,:),'linewidth',2)
                plot(f3_tot(:,plotIndices(ii)),G3(:,plotIndices(ii))*1e6,Marker3D,'Color',COMSOL3Colors(ii+COffset,:),'linewidth',C3LW)
                plot(f_tot(:,plotIndices(ii)),G2(:,plotIndices(ii))*1e6,Marker2D,'Color',COMSOLColors(ii+COffset,:),'linewidth',COMSOLLW(ii));
                
                %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
                LegStrings{4*(ii-1) +3} = sprintf('Model (Piston Impedance) %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{4*(ii-1) +2} = sprintf('COMSOL 3D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{4*(ii-1) +1} = sprintf('COMSOL 2D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{4*(ii-1) +4} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
                
            else
                %plot(test.f{plotIndices(ii)}/1e6,abs(test.xdef{plotIndices(ii)})*1e9,'-','Color',CustomColours(1+ii,:),'linewidth',1.75)
                plot(test.f{plotIndices(ii)}/1e6,real(test.G{plotIndices(ii)})*1e6,'-','Color',Model3Colors(ii+COffset,:),'linewidth',ModelLineW)
                plot(testP.f{plotIndices(ii)}/1e6,real(testP.G{plotIndices(ii)})*1e6,'-','Color',CustomColours(4,:),'linewidth',2)
                plot(f_tot(:,plotIndices(ii)),G2(:,plotIndices(ii))*1e6,Marker2D,'Color',COMSOLColors(ii+COffset,:),'linewidth',COMSOLLW(ii));
                plot(f3_tot(:,plotIndices(ii)),G3(:,plotIndices(ii))*1e6,Marker3D,'Color',COMSOL3Colors(ii+COffset,:),'linewidth',C3LW)

                %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
                LegStrings{4*(ii-1) +1} = sprintf('COMSOL 3D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{4*(ii-1) +2} = sprintf('COMSOL 2D %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{4*(ii-1) +4} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));
                LegStrings{4*(ii-1) +3} = sprintf('Model (Piston Impedance) %.0fV',test.Vop(plotIndices(ii)));
            end





        else
            %order swapped:
            plot(test.f{plotIndices(ii)}/1e6,real(test.G{plotIndices(ii)})*1e6,'-','Color',CustomColours(1+ii,:),'linewidth',2)
            plot(f_tot(:,plotIndices(ii)),G2(:,plotIndices(ii))*1e6,'-','Color',COMSOLColors(ii,:),'linewidth',COMSOLLW(ii));
            %these are swapped because they get swapped again later - it's a dumb solution - swap the order of the leg strings if not plotting in reverse order
            LegStrings{2*(ii-1) +1} = sprintf('COMSOL %.0fV',test.Vop(plotIndices(ii)));
            LegStrings{2*(ii-1) +2} = sprintf('Model %.0fV',test.Vop(plotIndices(ii)));

        end
    end
    %ylabel('Conductance (Ohms^-^1)')
    ylabel('Conductance (µΩ^-^1)')
    
    xlabel('Frequency (MHz)')
    
    %if reversed:
    LegStrings = LegStrings(length(LegStrings):-1:1);

    legend(LegStrings); 
    legend('boxoff');
    legend('location','southeast');

    ax1 = gca;
    ax1.FontSize = 14;
    ax1.LineWidth = 1.75;
    %ax1.Legend.FontSize = 14;
    ax1.FontName = 'Arial';
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    title('Model vs FEM (Oil)');
    
    if Zplot_LogScale
        set(gca, 'YScale', 'log')
    end 
    
    if exist("minYlim")
        if exist("maxYlim")
            ylim([minYlim maxYlim])
        else
            ylim([minYlim inf])
        end
    elseif exist("maxYlim")
        ylim([-inf maxYlim])
    end
    

    [~,maxIdxModel] = max(test.G{plotIndices(ii)});
    [~,maxIdxPiston] = max(testP.G{plotIndices(ii)});
    [~,maxIdx3D] = max(G3(:,plotIndices(ii)));
    [~,maxIdx2D] = max(G2(:,plotIndices(ii)));
    fmax_model = test.f{1}(maxIdxModel);
    fmax_piston = testP.f{1}(maxIdxPiston);
    fmax_3D = f3_tot(maxIdx3D,plotIndices(ii));
    fmax_2D = f_tot(maxIdx2D,plotIndices(ii));
    fprintf('\nFrequency of Peak Conductance\n')
    fprintf('Model: %.2f (MHz)\n',fmax_model/1e6);
    fprintf('Model (Piston): %.2f (MHz)\n',fmax_piston/1e6);
    fprintf('COMSOL 3D: %.2f (MHz)\n',fmax_3D);
    fprintf('COMSOL 2D: %.2f (MHz)\n',fmax_3D);  
end