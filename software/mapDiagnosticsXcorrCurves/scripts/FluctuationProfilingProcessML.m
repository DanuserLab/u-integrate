classdef FluctuationProfilingProcessML < DataProcessingProcessML
    % Process Class for Fluctuation Profiling Around Edge Motion Events (FPAEME)
    % FluctuationProfilingML.m is the wrapper function
    % FluctuationProfilingProcessML is part of XcorrFluctuationPackage, it requires MultiScaleAutoSegmentationProcess, MaskRefinementProcess,
    % ProtrusionProcess, WindowingProcess, ProtrusionSamplingProcess, and WindowSamplingProcess to be done.
    %
    % Changed this process name from "Fluctuation Profiling" to "Fluctuation Profiling Around Motion Events" on 10/23/2018.
    %
    % Qiongjing (Jenny) Zou, Sep 2018

    methods (Access = public)
        function obj = FluctuationProfilingProcessML(owner, varargin)
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
                super_args{2} = FluctuationProfilingProcessML.getName;
                super_args{3} = @FluctuationProfilingML;
                if isempty(funParams)
                    funParams=FluctuationProfilingProcessML.getDefaultParams(owner,outputDir);
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
            h= @FluctuationProfilingProcessMLGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieList'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            MDs = owner.getMovies();
            funParams.MovieDataIndex = 1:numel(MDs);
            % Listed all default parameters in phaseMasking, phaseDescriptives_OneChan, phaseDescriptives_MaxMinVel_OneChan, and MLsummary_FluctuationProfiling.
            
            % MD level params:
            % Used first MD's numel(MD.channels_) for ChannelIndex, however 
            % did Sanity Checks in FluctuationProfilingML.m to make sure numel(MD.channels_) consistent over all MDs
            load(owner.movieDataFile_{1}) % load sample MD.
            funParams.ChannelIndex = 1 : numel(MD.channels_);

            funParams.timeInterval = MD.timeInterval_; 

            % Did not register outputdir on MD level to this process:
            % funParams.OutputDirectory = [outputDir  filesep 'FluctuationProfiling'];

            funParams.smParamTh = 0.4; % default value was not given in phaseMasking.m, but in script.
            funParams.minimumRunLength = 5; % default is 1; 5 used in script.
            funParams.subFrames = [];
            funParams.impute = true;
            funParams.figFlag = 'off';
            % funParams.omittedWindows = [];
            funParams.omittedWindows = false; % created parameter when create GUI to replace original omittedWindows
            funParams.Folding = false;
            funParams.movingAvgSmoothing = true;

            funParams.maxLayer = 2;
            funParams.samplingBw = 20; % default value is not given in phaseDescriptives_OneChan.m, but in script.
            funParams.WithN = false;
             
            for i = 1:numel(MD.channels_)
                funParams.chanName{i} = ['channelName' num2str(i)]; % e.g. {'Actin' 'mDia1'}
                % funParams.OutputDirectory{i} = [outputDir  filesep 'Summary_PhaseDesc_Chan', num2str(i)]; % do not separate funParams.OutputDirectory here, otherwise folder icon on Control Panel of this package won't work.
            end
            funParams.OutputDirectory = outputDir; % ML level params.  % complete output folders see ML.processes_{x}.outFilePaths_

            funParams.lagMax0 = funParams.samplingBw; % 10/22/2018 Make = SamplingBw; if not nan, will use to calculate frSize in MLsummary_FluctuationProfiling.m
            
        end
    end
end