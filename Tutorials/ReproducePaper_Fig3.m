%Reproduce Fig 3 in our Paper:

%Add the path of the modeling code: (requires running the whole script)
filePath = matlab.desktop.editor.getActiveFilename; %get active file location
disp(filePath) %temp
TutorialsFolder = fileparts(filePath); ProjectFolder = fullfile(TutorialsFolder,'..');
CodePath = fullfile(ProjectFolder,"Code/");addpath(genpath(CodePath))

%define a struct with our device dimensions & analysis parameters
clear CMUT; CMUT = struct;
%device dimensions: 
CMUT.a = 60e-6; %membrane half-width [m]
CMUT.h = 5e-6; %membrane thickness [m]
CMUT.d0 = 420e-9; %cavity vacuum gap [m]
CMUT.b = 1.5e-3; %membrane half-length [m]
CMUT.di = 360e-9; %insulating layer thickness [m]

%material properties:
CMUT.epsr = 3.9; %relative permittivity of dielectric layer
CMUT.E = 148e9; %Young's Modulus of membrane [Pa]
CMUT.rho = 2329; %Density of membrane [kg/m^3]

%Metal electrode layer:
CMUT.h_metal = 0; %Electrode thickness [m]

test = edrc(CMUT); %simulate results

%% Processed FEM maximum deflection vs voltage:
FEMVW = zeros(16,2);
FEMVW(:,1) = [0:10:110 117 118 130 140];
FEMVW(:,2) = [
   34.2897
   34.9302
   36.8676
   40.1510
   44.8769
   51.1819
   59.2824
   69.5476
   82.4777
   99.1176
  121.2660
  155.2239
  213.3077
  420.1281
  420.1354
  420.1190];

%% Plot settings:

%Figure settings
ScreenSize = get(0,'ScreenSize');
FigWidth = 500;
FigHeight = 400;
XStartPos = 100;
YStartPos = 595;
NoFigs = 0;
figoffset = 100;

CustomColors = 1;
CustomColours(1,:) = 0.15*[1 1 1];
CustomColours(2,:) = 1/255 * [72 118 255]; %blue
CustomColours(3,:) = 1/255 * [220 20 60]; %red
CustomColours(4,:) = 1/255 * [255 215 0]*0.9; %Gold
CustomColours(5,:) = 1/255 * [139 0 139]; %purple
CustomColours(6,:) = 1/255 * [0 205 102]; %green

%plot the figure:
VWFig = figure(figoffset+NoFigs);
clf;
try
    VWFig.Position = GenerateFigPositions(NoFigs,FigWidth,FigHeight, XStartPos,0,ScreenSize,85);
catch
    VWFig.Position = [XStartPos YStartPos FigWidth FigHeight];
end
NoFigs = NoFigs + 1;
hold on
box on

plot(test.VDC,test.w0*1e9,'Color',CustomColours(1,:),'linewidth',2);
xlabel('Bias Voltage (V)');
ylabel('Peak deflection w_0 (nm)');
title('Membrane Deflection with Increasing Voltage');
plot(FEMVW(:,1),FEMVW(:,2),':sq','Color', CustomColours(2,:),'linewidth',2);
legend('1D Model','FEM');
legend('location','northwest');
legend('boxoff')

%for shrinking fig:
ax1 = gca;
ax1.FontSize = 14;

ax1.LineWidth = 1.75;
ax1.Legend.FontSize = 14;
ax1.FontName = 'Arial';

ax1.XColor = 'k';
ax1.YColor = 'k';



