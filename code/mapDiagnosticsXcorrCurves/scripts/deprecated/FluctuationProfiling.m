function FluctuationProfiling(movieDataOrProcess, varargin)
% FluctuationProfiling wrapper function for phaseMasking, phaseDescriptives_OneChan,
% phaseDescriptives_MaxMinVel_OneChan, (and MLsummary_FluctuationProfiling). to be executed by
% FluctuationProfilingProcess.
% NOTE: MLsummary_FluctuationProfiling is not executed here, since it reqires ML as input.
%
% INPUT
% movieDataOrProcess - either a MovieData (legacy)
%                      or a Process (new as of July 2016)
%
% param - (optional) A struct describing the parameters, overrides the
%                    parameters stored in the process (as of Aug 2016)
%
% OUTPUT
% none (saved to p.OutputDirectory)
%
% Changes
% As of July 2016, the first argument could also be a Process. Use
% getOwnerAndProcess to simplify compatability.
%
% As of August 2016, the standard second argument should be the parameter
% structure
%
% Qiongjing (Jenny) Zou, Aug 2018


%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('param',[], @isstruct);
ip.parse(movieDataOrProcess, varargin{:});
p = ip.Results;
paramsIn = p.param;

%% Registration

% Get MovieData object and Process
% If movieDataOrProcess is a MovieData and does not contain an
% FluctuationProfilingProcess, create FluctuationProfilingProcess using constructor with no
% arguments.
% If movieDataOrProcess is a MovieData and does contain an FluctuationProfilingProcess,
% then return the first instance of an FluctuationProfilingProcess.
% If movieDataOrProcess is an FluctuationProfilingProcess, then return the Process and it's
% MovieData owner.
% Otherwise throw an error.
[movieData, process] = getOwnerAndProcess(movieDataOrProcess,'FluctuationProfilingProcess',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used 
% rather than the one stored in FluctuationProfilingProcess

% Parameters
currSmParamTh = p.smParamTh;
currMinimumRunLength = p.minimumRunLength;
currSubFrames = p.subFrames;
currImpute = p.impute;
currFigFlag = p.figFlag;
currFolding = p.Folding;
currMovingAvgSmoothing = p.movingAvgSmoothing;
currMaxLayer = p.maxLayer;
currChanName = p.chanName; % cell array, # cell = # channels
currSamplingBw = p.samplingBw;
currWithN = p.WithN;

% currLagMax = p.lagMax;
% currTimeInterval = p.timeInterval;

    % If LBTestProcess is done, exclude quiescent windows (calculate currMovingAvgSmoothing).
LBTestId = movieData.getProcessIndex('LBTestProcess');
if ~isempty(LBTestId)
    LBTestOutDir = movieData.processes_{LBTestId}.outFilePaths_;
    indPath = fullfile(LBTestOutDir, 'indActive_windowIndex.mat');
    raw = load(indPath);

    omitWin = find(raw.indActive == 0);
    disp('== omittedWindows: ==')
    disp(omitWin)
    currOmittedWindows = omitWin;
else
    currOmittedWindows = p.omittedWindows;
end


% Sanity Checks
nChan = numel(movieData.channels_);
if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

% precondition / error checking:
windowingId = movieData.getProcessIndex('WindowingProcess');
winProc = movieData.getProcess(windowingId);
methodType = winProc.funParams_.MethodName;
if ~isequal(methodType, 'ConstantNumber')
    error("Method used to propagate the windows from one frame to the next needs to be ConstantNumber")
end

protSamplingId = movieData.getProcessIndex('ProtrusionSamplingProcess');
if isempty(protSamplingId)
    error("ProtrusionSamplingProcess needs to be done before run this process.")
end

winSamplingId = movieData.getProcessIndex('WindowSamplingProcess');
if isempty(winSamplingId)
    error("WindowSamplingProcess needs to be done before run this process.")
end

if isempty(movieData.pixelSize_); error('movieData.pixelSize_ is required.'); end
if isempty(movieData.timeInterval_); error('movieData.timeInterval_ is required.'); end

% logging input paths (bookkeeping)
nChan = numel(movieData.channels_);
inFilePaths = cell(3, nChan);
for i = p.ChannelIndex
    protSamplingProc = movieData.getProcess(protSamplingId);
    inFilePaths{1,i} = protSamplingProc.outFilePaths_;% there is only one outFilePaths_ in ProtrusionSamplingProcess for all channels.

    winSamplingProc = movieData.getProcess(winSamplingId);
    inFilePaths{2,i} = winSamplingProc.outFilePaths_{1,i};

    LBTestId = movieData.getProcessIndex('LBTestProcess');
    if ~isempty(LBTestId)
        LBTestProc = movieData.getProcess(LBTestId);
        inFilePaths{3,i} = LBTestProc.outFilePaths_;% there is only one outFilePaths_ in LBTestProcess for all channels.
    else 
        inFilePaths{3,i} = [];
    end
end
process.setInFilePaths(inFilePaths);

% logging output paths
% Did not separate MD level output for channels.
mkClrDir(p.OutputDirectory);
process.setOutFilePaths(p.OutputDirectory);

%% Algorithm

for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    chanName = currChanName{i};

    % phaseMasking:
    phaseDescDirName = 'phaseDescriptives';
    figuresDir = fullfile(p.OutputDirectory, phaseDescDirName);

    [protMask1, retMask1] = phaseMasking(movieData, currSmParamTh, figuresDir, 'minimumRunLength', currMinimumRunLength, ...
            'subFrames', currSubFrames,'impute', currImpute, 'figFlag', currFigFlag, 'omittedWindows', currOmittedWindows, ...
            'Folding', currFolding, 'movingAvgSmoothing', currMovingAvgSmoothing);

    % phaseDescriptives_OneChan:
        % Retraction (use retMask1):
    chanNameTag = [chanName, '-ret'];
    chanTitle = chanNameTag;

    phaseDescriptives_OneChan(movieData, iChan, currMaxLayer, chanNameTag, chanTitle, retMask1, currSamplingBw, figuresDir, ...
            'WithN', currWithN, 'impute', currImpute, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
            'movingAvgSmoothing', currMovingAvgSmoothing, 'figFlag', currFigFlag, 'Folding', currFolding)

        % Protrusion (use protMask1):
    chanNameTag = [chanName, '-prot'];
    chanTitle = chanNameTag;

    phaseDescriptives_OneChan(movieData, iChan, currMaxLayer, chanNameTag, chanTitle, protMask1, currSamplingBw, figuresDir, ...
            'WithN', currWithN, 'impute', currImpute, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
            'movingAvgSmoothing', currMovingAvgSmoothing, 'figFlag', currFigFlag, 'Folding', currFolding)

    % phaseDescriptives_MaxMinVel_OneChan:
    chanTitle = chanName;

    phaseDescriptives_MaxMinVel_OneChan(movieData, iChan, currMaxLayer, chanName, chanTitle, currSmParamTh, currSamplingBw, ...
            figuresDir, 'WithN', currWithN, 'impute', currImpute, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
            'movingAvgSmoothing', currMovingAvgSmoothing, 'figFlag', currFigFlag, 'Folding', currFolding)
        
    close all

    % MLsummary_FluctuationProfiling
    % summaryDirName = 'SummaryFPAEME';

    % Continue here, think about how to accommodate ML into the process
    % MLsummary_FluctuationProfiling(ML, chanName, currMaxLayer, phaseDescDirName, ...
    %     summaryDirName, 'lagMax0', currLagMax, 'timeInterval', currTimeInterval)
    % close all

end

end