classdef FluctuationProfilingProcess < DataProcessingProcess
    % Process Class for Fluctuation Profiling Around Edge Motion Events (FPAEME)
    % FluctuationProfiling.m is the wrapper function
    % FluctuationProfilingProcess is part of XcorrFluctuationPackage, it requires MultiScaleAutoSegmentationProcess, MaskRefinementProcess,
    % ProtrusionProcess, WindowingProcess, ProtrusionSamplingProcess, and WindowSamplingProcess to be done.
    % TODO: now write for MD, will need to think a way to accommodate ML.
    %
    % Changed this process name from "Fluctuation Profiling" to "Fluctuation Profiling Around Motion Events" on 10/23/2018.
    %
    % Qiongjing (Jenny) Zou, Aug 2018

    methods (Access = public)
        function obj = FluctuationProfilingProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else 
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = FluctuationProfilingProcess.getName;
                super_args{3} = @FluctuationProfiling;
                if isempty(funParams)
                    funParams=FluctuationProfilingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});

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
            h= @cliGUI; % TODO Generic Command Line Interface GUI, need to change to a customized FluctuationProfilingProcessGUI.
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % Listed all default parameters in phaseMasking, phaseDescriptives_OneChan, phaseDescriptives_MaxMinVel_OneChan, (and MLsummary_FluctuationProfiling).
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'FluctuationProfiling'];

            funParams.smParamTh = 0.4 % default value was not given in phaseMasking.m, but in script.
            funParams.minimumRunLength = 5 % default is 1; 5 used in script.
            funParams.subFrames = [];
            funParams.impute = true;
            funParams.figFlag = 'off';
            funParams.omittedWindows = [];
            funParams.Folding = false;
            funParams.movingAvgSmoothing = false;

            funParams.maxLayer = 2;
            funParams.samplingBw = 20; % default value is not given in phaseDescriptives_OneChan.m, but in script.
            funParams.WithN = false;
             
            for i = 1:numel(owner.channels_)
                funParams.chanName{i} = ['channelName' num2str(i)]; % e.g. 'Actin' or 'mDia1'
            end

            % funParams.lagMax0 = 5; % fake value, if lagMax0 is default (nan) or other specified, will calculated in MLsummary_FluctuationProfiling.m;
            % funParams.timeInterval = MovieData.timeInterval_;
            
        end
    end
end