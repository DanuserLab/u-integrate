classdef FluctuationProfilingForGCAProcessML < DataProcessingProcessML
    % Process Class for Fluctuation Profiling Around Edge Motion Events (FPAEME)
    %
    % Updates:
    % J Noh, 7/2021. Modified for GCA package.
    % Qiongjing (Jenny) Zou, Sep 2018
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
%
% This file is part of GrangerCausalityAnalysisPackage.
% 
% GrangerCausalityAnalysisPackage is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GrangerCausalityAnalysisPackage is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GrangerCausalityAnalysisPackage.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

    methods (Access = public)
        function obj = FluctuationProfilingForGCAProcessML(owner, varargin)
            if nargin == 0
                super_args = {};
            else 
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar); % changed to ML's outputDir
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = FluctuationProfilingForGCAProcessML.getName;
                super_args{3} = @FluctuationProfilingForGCAMLWrapper;
                if isempty(funParams)
                    funParams=FluctuationProfilingForGCAProcessML.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcessML(super_args{:});

        end
        
        
        % function output = loadChannelOutput(obj, iChan, varargin)
        %     outputList = {};
        %     nOutput = length(outputList);

        %     ip.addRequired('iChan',@(x) obj.checkChanNum(x));
        %     ip.addOptional('iOutput',1,@(x) ismember(x,1:nOutput));
        %     ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
        %     ip.addParamValue('useCache',false,@islogical);
        %     ip.parse(iChan,varargin{:})
    
        %     s = cached.load(obj.outFilePaths_{iChan},'-useCache',ip.Results.useCache);

        %     output = s.Imean;          
        % end
    end

    methods (Static)
        function name = getName()
            name = 'Fluctuation Profiling Around Motion Events';
        end

        function h = GUI()
            h= @FluctuationProfilingForGCAProcessMLGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieList'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set up global parameters from SetUpGCAParametersProcessML
            setupProcID = owner.getProcessIndex('SetUpGCAParametersProcessML');
            funParams = parseProcessParams(owner.processes_{setupProcID}); 
            
            % set the frame size for moving median normalization
            if funParams.lowFreqSubtraction      
                funParams.movMedFrameSize = funParams.avgLengthOfProtRetCycle;
            else
                funParams.movMedFrameSize = [];
            end 
             
            % spline smoothing param (smParam) for Xcorr map is not used.
            funParams.excludeQuiescentWindows = true; 
            
            % For FPAME
            funParams.smParamPhaseSegmentation = 0.5; % smParam used to segment Prot/Ret phases
            funParams.minimumRunLength = 3; % minimal frame size for P or R phases
            % +-samplingBw frames are sampled around an edge motion event
            % used for phaseDescriptives_OneChan.m and
            % MLsummary_FluctuationProfiling.m
            funParams.samplingBw = funParams.avgLengthOfProtRetCycle; 
            
        end
    end
end