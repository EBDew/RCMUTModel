function [Cout] = edrc(InStruct)
    %This function is used to call the RCMUTModel() function with an structure of parameters.
    %If you try to call RCMUTModel(struct), the argument parser will throw an error because the struct is not a parameter.
    %This function converts an input struct into a set of name value pairs to be used with RCMUTModel, then returns the output.
    %This is the most convenient way to call RCMUTModel.
    %For a list of supported parameters, and their default values either visit the documentation for RCMUTModel, or type the following into the console:
    %>>test = RCMUTModel();test.Params %this runs the code with its default parameters, then returns a structure with all of the parameters the code ran with
    %This function also includes some code to remove legacy parameters that are no longer supported, and convert them to the current parameter name
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

    % Allow this function to run with no or empty arguments
    if nargin == 0 || isempty(InStruct)
        Cout = RCMUTModel();
        return
    end
   
    
    %Convert some legacy parameters that existed in previous implementations of this function for backwards compatibility.

    %% Remove Parameters for Deprecated Features
    
    if isfield(InStruct,'CMVolumeFraction') 
        InStruct = rmfield(InStruct,'CMVolumeFraction'); %this parameter is deprecated, just delete it.
    end

    if isfield(InStruct,'scaleZbyAlpha') 
        InStruct = rmfield(InStruct,'scaleZbyAlpha'); %this parameter is deprecated, just delete it.
    end

    if isfield(InStruct,'scaleCompliance') 
        InStruct = rmfield(InStruct,'scaleCompliance'); %this parameter is deprecated, just delete it.
    end

    if isfield(InStruct,'h2')
        InStruct = rmfield(InStruct,'h2'); %this parameter is not needed, just delete it.
    end

    if isfield(InStruct,'calcImpedance')
        InStruct = rmfield(InStruct,'calcImpedance'); %deprecated
    end

    if isfield(InStruct,'LoadRadiationImpedance')
        InStruct = rmfield(InStruct,'LoadRadiationImpedance'); %deprecated
    end
   

    if isfield(InStruct,'SymbolicPostPrecision') 
        InStruct = rmfield(InStruct,'SymbolicPostPrecision'); %this parameter is deprecated, just delete it.
    end

    if isfield(InStruct,'calcDynamic') %deprecated parameter
        InStruct = rmfield(InStruct,'calcDynamic');
    end

    if isfield(InStruct,'calcDynamicModelParameters') %deprecated parameter
        InStruct = rmfield(InStruct,'calcDynamicModelParameters');
    end

    %% Backwards Compatibility: Deal with parameters whose names have been changed

    if isfield(InStruct,'DefaultFsweep') %this Parameter used to have a less descriptive name:
        InStruct.UseDefaultFrequencySweep = InStruct.DefaultFsweep;
        InStruct = rmfield(InStruct,'DefaultFsweep');
    end

    if isfield(InStruct,'CustomSweep') %this Parameter used to have a less descriptive name:
        InStruct.CustomFrequencySweepRange = InStruct.CustomSweep;
        InStruct = rmfield(InStruct,'CustomSweep');
    end

    if isfield(InStruct,'fstepP') %this Parameter used to have a less descriptive name: (precise frequency step)
        InStruct.fstepFine= InStruct.fstepP;
        InStruct = rmfield(InStruct,'fstepP');
    end

    %more descriptive names for circuit element parameters:
    if isfield(InStruct,'SeriesCapacitor') %this Parameter used to have a less descriptive name:
        InStruct.EnableSeriesCapacitor = InStruct.SeriesCapacitor;
        InStruct = rmfield(InStruct,'SeriesCapacitor');
    end

    %more descriptive names for circuit element parameters:
    if isfield(InStruct,'SeriesResistor') %this Parameter used to have a less descriptive name:
        InStruct.EnableSeriesResistor = InStruct.SeriesResistor;
        InStruct = rmfield(InStruct,'SeriesResistor');
    end

    %more descriptive names for circuit element parameters:
    if isfield(InStruct,'ParasiticCap') %this Parameter used to have a less descriptive name:
        InStruct.EnableParasiticCapacitance = InStruct.ParasiticCap;
        InStruct = rmfield(InStruct,'ParasiticCap');
    end
    
    %more descriptive names:
    if isfield(InStruct,'CircularPosts') %this Parameter used to have a less descriptive name:
        InStruct.circularPostAreaCorrection  = InStruct.CircularPosts;
        InStruct = rmfield(InStruct,'CircularPosts');
    end

    %more descriptive names:
    if isfield(InStruct,'PistonZ') %this Parameter used to have a less descriptive name:
        InStruct.UsePistonRadiationImpedance = InStruct.PistonZ;
        InStruct = rmfield(InStruct,'PistonZ');
    end

    %more descriptive names:
    if isfield(InStruct,'calcHyst') %this Parameter used to have a less descriptive name:
        InStruct.calcHysteresis = InStruct.calcHyst;
        InStruct = rmfield(InStruct,'calcHyst');
    end

    %more descriptive names: Radiation Impedance Source:
    if isfield(InStruct,'ClampedRadZSource') %this Parameter used to have a less descriptive name:
        switch InStruct.ClampedRadZSource
            case {2,3} %this used to refer to the 1D case
                InStruct.Clamped1DRadiationImpedance = 1;
            case {1,0} %this used to refer to the 2D case
                InStruct.Clamped1DRadiationImpedance = 2;
            otherwise %assume 1D
                InStruct.Clamped1DRadiationImpedance = 1;
        end
        InStruct = rmfield(InStruct,'ClampedRadZSource');
    end

    if isfield(InStruct,'a2')
        InStruct.a = InStruct.a2; %a was called 'a2' in some old versions of this code
        InStruct = rmfield(InStruct,'a2');
    end

    if isfield(InStruct,'h1')
        InStruct.h = InStruct.h1; %h was called 'h1' in some old versions of this code
        InStruct = rmfield(InStruct,'h1');
    end

    %% Parameters with Abbreviated names:
    %these can be used to save some typing
    if isfield(InStruct,'ZradD')
                InStruct.Clamped1DRadiationImpedance = InStruct.ZradD;
        InStruct = rmfield(InStruct,'ZradD');
    end

    if isfield(InStruct,'atr')
                InStruct.IsolatingTrenchWidth = InStruct.atr;
        InStruct = rmfield(InStruct,'atr');
    end

    
%% Call the function

    %convert the input structure to name value pairs and call function:
    nvp = namedargs2cell(InStruct); %convert to a name value pair of the parameters and their values
    Cout = RCMUTModel(nvp{:}); %old version
end
