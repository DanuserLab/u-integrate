classdef XcorrAnalysisProcess < DataProcessingProcess
    % Process Class for map Descriptive, Diagnotics and Xcorrelation Analysis.
    % XcorrAnalysis.m is the wrapper function
    % XcorrAnalysisProcess is part of XcorrFluctuationPackage, it requires MultiScaleAutoSegmentationProcess, MaskRefinementProcess,
    % ProtrusionProcess, WindowingProcess, ProtrusionSamplingProcess, and WindowSamplingProcess to be done.
    % TODO: now write for MD, will need to think a way to accommodate ML.
    %
    % Qiongjing (Jenny) Zou, Aug 2018

    methods (Access = public)
        function obj = XcorrAnalysisProcess(owner, varargin)
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
                super_args{2} = XcorrAnalysisProcess.getName;
                super_args{3} = @XcorrAnalysis;
                if isempty(funParams)
                    funParams=XcorrAnalysisProcess.getDefaultParams(owner,outputDir);
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
            name = 'Cross Correlation Analysis';
        end

        function h = GUI()
            h= @cliGUI; % TODO Generic Command Line Interface GUI, need to change to a customized XcorrAnalysisProcessGUI.
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % Listed all default parameters in mapDescriptives_OneChan, mapXcorrCurvePermutation_Vel, (and MLsummary_XcorrCurvesVelAcf)
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'Xcorr_Analysis'];

            % missing ichan 0 or 1
            funParams.maxLayer = 2; % for ichan=0, needs to be 1.
            funParams.adf = true; % default is false
            funParams.impute = false;
            funParams.movingAvgSmoothing = false;
            funParams.omittedWindows = [];
            funParams.subFrames = [];
            funParams.numPerm = 100; % default is 1000
            funParams.topograph = 'on';
            funParams.figFlag = 'off';
            funParams.WithN = false;
            funParams.parpoolNum = 4;
            funParams.rseed = 'shuffle';
            funParams.Folding = false;
            funParams.smParam = 0.8;

            funParams.lagMax = 5; % fake value, if lagMax is default (not specified), will calculated using round(tmax/4) in mapXcorrCurvePermutation_Vel.m or mapXcorrCurvePermutation.m;
            funParams.fullRange = false;
            funParams.h0 = [];
            funParams.mvFrSize = 0;

            for i = 1:numel(owner.channels_)
                funParams.chanName{i} = ['channelName' num2str(i)]; % e.g. 'Actin' or 'mDia1'
            end

%             funParams.timeInterval = MovieData.timeInterval_;
        end
    end
end