classdef QuiescentWindowDetectionProcessML < DataProcessingProcessML
    % Process Class for identifying quiescent windows by using Ljung-Box test.
    % QuiescentWindowDetectionML.m is the wrapper function
    % QuiescentWindowDetectionProcessML is part of XcorrFluctuationPackage, it requires MultiScaleAutoSegmentationProcess, MaskRefinementProcess,
    % ProtrusionProcess, WindowingProcess, ProtrusionSamplingProcess, and WindowSamplingProcess to be done.
    % 
    % Changed this process name from "LB Test" to "Quiescent Window Detection" on 10/23/2018.
    %
    % QuiescentWindowDetectionProcessML.m was created to replace LBTestProcessML.m on 10/24/2018.
    %
    % Qiongjing (Jenny) Zou, Aug 2018

    methods (Access = public)
        function obj = QuiescentWindowDetectionProcessML(owner, varargin)
            if nargin == 0
                super_args = {};
            else 
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);  % changed to ML's outputDir
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = QuiescentWindowDetectionProcessML.getName;
                super_args{3} = @QuiescentWindowDetectionML;
                if isempty(funParams)
                    funParams=QuiescentWindowDetectionProcessML.getDefaultParams(owner,outputDir);
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
            name = 'Quiescent Window Detection';
        end

        function h = GUI()
            h= @QuiescentWindowDetectionProcessMLGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieList'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % Listed all default parameters in mapDescriptives_Vel_LB.
            % funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.ProtSamplingProcessIndecx = []; % specify which protsamplingProcess to use, default is the last one.

            funParams.impute = true;
            funParams.movingAvgSmoothing = true;
            funParams.figFlag = 'off';
            funParams.omittedWindows = [];
            funParams.subFrames = [];
            funParams.topograph = 'on';
            funParams.Folding = false;
            funParams.derivative = false;
            funParams.smParam = 0.8;
            funParams.movmeanNum = 3;

            funParams.ratioThreshold = 0.5;

            funParams.OutputDirectory = outputDir; % ML level params.  % complete output folders see ML.processes_{x}.outFilePaths_
        end
    end
end