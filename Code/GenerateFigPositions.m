function [FigPosition] = GenerateFigPositions(NumFigs,FigWidth,FigHeight, Xoffset,Yoffset,ScreenSize,FigToolBarHeight)
% This function is used to generate the positions of a grid of equally-sized figures that fits on a given monitor.
% NumFigs - Number of figures already plotted (starts at 0)
% FigWidth - Width of each figure (in pixels)
% FigHeight - Height of each figure (in pixels)
% Xoffset - Offset from left side of the screen for the leftmost figure (in pixels)
% Yoffset - Offset from top of the screen for the top figure (in pixels)
% ScreenSize - vector of [1 1 NumberOfXPixels NumberOfYPixels] (in pixels). eg [1 1 1920 1080]. This is optional as this function can calculate it, however,
%   entering it as an arugment may improve performance
% FigToolBarHeight - Height of the MATLAB toolbar above each figure in (in pixels).
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
    NumFigs (1,1) {mustBeInteger,mustBeNonnegative} = 0 %number of figures already plotted
    FigWidth (1,1) {mustBePositive} = 500 %Width of the figure (in pixels)
    FigHeight (1,1) {mustBePositive} = 400 %Height of the figure (in pixels)
    Xoffset (1,1) {mustBeNonnegative} = 100 %Offset from left side of the screen (in pixels)
    Yoffset (1,1) {mustBeNonnegative} = 0 %Offset from top of the screen (in pixels)
    ScreenSize (1,4) {mustBePositive} = get(0,'ScreenSize')
    FigToolBarHeight (1,1) {mustBePositive} = 85 %Height of the toolbar above the figure (in pixels)
end
AvailableHeight = ScreenSize(4) - Yoffset;
TotalFigHeight = FigHeight + FigToolBarHeight;
%determine how many figures fit along each dimension
NumberOfXFigs = floor( (ScreenSize(3)-Xoffset)/ FigWidth );
NumberOfYFigs = floor( AvailableHeight/ TotalFigHeight);
MaximumNumberOfFigs = NumberOfXFigs*NumberOfYFigs;

%loop back to the start if there are more figures than fit on a monitor:
%Eg for a monitor that fits 3 figures wide and 2 tall, fig 7 is in the same place as fig 1
if NumFigs >= MaximumNumberOfFigs
    NumFigs = mod(NumFigs,MaximumNumberOfFigs);
end

YStartPosition = AvailableHeight - (floor(NumFigs/NumberOfXFigs)+1)* TotalFigHeight; %Position of the bottom of the first fig.
XStartPosition = Xoffset + mod(NumFigs,NumberOfXFigs)*FigWidth; %Position of the left edge of the first fig.

FigPosition = [XStartPosition YStartPosition FigWidth FigHeight];
end
