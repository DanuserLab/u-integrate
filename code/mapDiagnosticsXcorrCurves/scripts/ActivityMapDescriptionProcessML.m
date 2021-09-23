classdef ActivityMapDescriptionProcessML < DataProcessingProcessML
    % Process Class for map Descriptive.
    % ActivityMapDescriptionML.m is the wrapper function
    % ActivityMapDescriptionProcessML is part of XcorrFluctuationPackage, it requires MultiScaleAutoSegmentationProcess, MaskRefinementProcess,
    % ProtrusionProcess, WindowingProcess, ProtrusionSamplingProcess, and WindowSamplingProcess to be done.
    %
    % On 10/30/2018 ActivityMapDescriptionProcessML was separated from XcorrAnalysisProcessML as
    % per Jungsik's request.
    %
    % Qiongjing (Jenny) Zou, Oct 2018

    methods (Access = public)
        function obj = ActivityMapDescriptionProcessML(owner, varargin)
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
                super_args{2} = ActivityMapDescriptionProcessML.getName;
                super_args{3} = @ActivityMapDescriptionML;
                if isempty(funParams)
                    funParams=ActivityMapDescriptionProcessML.getDefaultParams(owner,outputDir);
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
            name = 'Activity Map Description';
        end

        function h = GUI()
            h= @ActivityMapDescriptionProcessMLGUI;
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
            % Listed all default parameters in mapDescriptives_OneChan.

            % Used first MD's numel(MD.channels_) for ChannelIndex, however 
            % did Sanity Checks in ActivityMapDescriptionML.m to make sure numel(MD.channels_) consistent over all MDs
            load(owner.movieDataFile_{1}) % load sample MD.
            funParams.ChannelIndex = 1 : numel(MD.channels_);
            
            funParams.OutputDirectory = outputDir; % ML level params.  % complete output folders see ML.processes_{x}.outFilePaths_

            funParams.maxLayer = 2; % for ichan=0, needs to be 1.
            funParams.adf = true; % default is false
            funParams.impute = true;
            funParams.movingAvgSmoothing = true;
            % funParams.omittedWindows = [];
            funParams.omittedWindows = false; % created parameter when create GUI to replace original omittedWindows
            funParams.subFrames = [];
            funParams.numPerm = 100; % default is 1000
            funParams.topograph = 'on';
            funParams.figFlag = 'off';
            funParams.WithN = false;
            % funParams.parpoolNum = 4; % parpoolNum is no longer a needed parameter for mapDescriptives_OneChan
            funParams.rseed = 'shuffle'; % no other options for this rseed parameter.
            funParams.Folding = false;
            funParams.smParam = 0.8;

            for i = 1:numel(MD.channels_)
                funParams.chanName{i} = ['channelName' num2str(i)]; % e.g. {'Actin' 'mDia1'}
            end
        end
    end
end