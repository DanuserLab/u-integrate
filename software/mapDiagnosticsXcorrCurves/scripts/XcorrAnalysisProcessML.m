classdef XcorrAnalysisProcessML < DataProcessingProcessML
    % Process Class for Xcorrelation Analysis.
    % XcorrAnalysisML.m is the wrapper function
    % XcorrAnalysisProcessML is part of XcorrFluctuationPackage, it requires MultiScaleAutoSegmentationProcess, MaskRefinementProcess,
    % ProtrusionProcess, WindowingProcess, ProtrusionSamplingProcess, and WindowSamplingProcess to be done.
    %
    % On 10/30/2018 ActivityMapDescriptionProcessML was separated from XcorrAnalysisProcessML as
    % per Jungsik's request.
    %
    % Qiongjing (Jenny) Zou, Aug 2018

    methods (Access = public)
        function obj = XcorrAnalysisProcessML(owner, varargin)
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
                super_args{2} = XcorrAnalysisProcessML.getName;
                super_args{3} = @XcorrAnalysisML;
                if isempty(funParams)
                    funParams=XcorrAnalysisProcessML.getDefaultParams(owner,outputDir);
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
            name = 'Cross Correlation Analysis';
        end

        function h = GUI()
            h= @XcorrAnalysisProcessMLGUI;
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
            % Listed all default parameters in mapDescriptives_OneChan, mapXcorrCurvePermutation_Vel, and MLsummary_XcorrCurvesVelAcf

            % MD level params:
            % Used first MD's numel(MD.channels_) for ChannelIndex, however 
            % did Sanity Checks in XcorrAnalysisML.m to make sure numel(MD.channels_) consistent over all MDs
            load(owner.movieDataFile_{1}) % load sample MD.
            funParams.ChannelIndex = 1 : numel(MD.channels_);

            funParams.timeInterval = MD.timeInterval_; 

            % Did not register outputdir on MD level to this process:
            % mapDescriptivesDirName = 'mapDescriptives';
            % funParams.OutputDirectory{1} = [outputDir  filesep 'Xcorr_Analysis' filesep mapDescriptivesDirName];

            % ML level params:
            funParams.OutputDirectory = outputDir; % do not separate funParams.OutputDirectory here, otherwise folder icon on Control Panel of this package won't work.
                                                % complete output folders see ML.processes_{x}.outFilePaths_
            % v = 1: numel(MD.channels_)+1; % +1 to include ch0 
            % combi = nchoosek(v,2);
            % [m, ~] = size(combi); % m possible combinations
            % for i = 1:m
            %     funParams.OutputDirectory{i} = [outputDir  filesep 'Summary_Xcorr_ch' num2str(combi(i,2)-1) 'ch' num2str(combi(i,1)-1)];
            % end

            funParams.maxLayer = 2; % for ichan=0, needs to be 1.
            funParams.impute = true;
            funParams.movingAvgSmoothing = true;
            % funParams.omittedWindows = [];
            funParams.omittedWindows = false; % created parameter when create GUI to replace original omittedWindows
            funParams.subFrames = [];
            funParams.numPerm = 100; % default is 1000
            funParams.topograph = 'on';
            funParams.figFlag = 'off';
            funParams.WithN = false;
            % funParams.parpoolNum = 4; % parpoolNum is no longer a needed parameter for mapXcorrCurvePermutation_Vel
            funParams.rseed = 'shuffle'; % no other options for this rseed parameter.
            funParams.Folding = false;

            funParams.lagMax = 5; % fake value, if lagMax is default (not specified), will calculated using round(tmax/4) in mapXcorrCurvePermutation_Vel.m or mapXcorrCurvePermutation.m;
            funParams.fullRange = false;
            funParams.h0 = []; % This parameter is not used in the function(s).
            funParams.mvFrSize = 0;

            for i = 1:numel(MD.channels_)
                funParams.chanName{i} = ['channelName' num2str(i)]; % e.g. {'Actin' 'mDia1'}
            end

        end
    end
end