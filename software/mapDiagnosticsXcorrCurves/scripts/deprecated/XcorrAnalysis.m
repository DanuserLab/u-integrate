function XcorrAnalysis(movieDataOrProcess, varargin)
% XcorrAnalysis wrapper function for mapDescriptives_OneChan, 
% mapXcorrCurvePermutation_Vel, (and MLsummary_XcorrCurvesVelAcf) to be executed by
% XcorrAnalysisProcess.
% NOTE: MLsummary_XcorrCurvesVelAcf is not executed here, since it reqires ML as input.
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
% XcorrAnalysisProcess, create XcorrAnalysisProcess using constructor with no
% arguments.
% If movieDataOrProcess is a MovieData and does contain an XcorrAnalysisProcess,
% then return the first instance of an XcorrAnalysisProcess.
% If movieDataOrProcess is an XcorrAnalysisProcess, then return the Process and it's
% MovieData owner.
% Otherwise throw an error.
[movieData, process] = getOwnerAndProcess(movieDataOrProcess,'XcorrAnalysisProcess',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used 
% rather than the one stored in XcorrAnalysisProcess

% Parameters
currMaxLayer = p.maxLayer;
currAdf = p.adf;
currImpute = p.impute;
currMovingAvgSmoothing = p.movingAvgSmoothing;
currSubFrames = p.subFrames;
currNumPerm = p.numPerm;
currTopograph = p.topograph;
currFigFlag = p.figFlag;
currWithN = p.WithN;
currParpoolNum = p.parpoolNum;
currRseed = p.rseed;
currFolding = p.Folding;
currSmParam = p.smParam;
currLagMax = p.lagMax;
currFullRange = p.fullRange;
currH0 = p.h0;
currMvFrSize = p.mvFrSize;
currChanName = p.chanName; % cell array, # cell = # channels
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

% mapDescriptives_OneChan:
mapDescriptivesDirName = 'mapDescriptives';
figuresDir1 = fullfile(p.OutputDirectory, mapDescriptivesDirName);

iChan0 = 0;
maxLayer = 1;
chan0Name = 'Vel';
chan0Title = 'Velocity (nm/sec)';
mapDescriptives_OneChan(movieData, iChan0, maxLayer, chan0Name, chan0Title, figuresDir1, 'adf', currAdf,'impute', currImpute, ...
    'movingAvgSmoothing', currMovingAvgSmoothing, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
    'numPerm', currNumPerm, 'topograph', currTopograph, 'figFlag', currFigFlag, 'WithN', currWithN, 'parpoolNum', currParpoolNum, ...
    'rseed', currRseed, 'Folding', currFolding, 'smParam', currSmParam)


for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    chanName = currChanName{i};
    chanTitle = chanName;

    mapDescriptives_OneChan(movieData, iChan, currMaxLayer, chanName, chanTitle, figuresDir1, 'adf', currAdf,'impute', currImpute, ...
        'movingAvgSmoothing', currMovingAvgSmoothing, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
        'numPerm', currNumPerm, 'topograph', currTopograph, 'figFlag', currFigFlag, 'WithN', currWithN, 'parpoolNum', currParpoolNum, ...
        'rseed', currRseed, 'Folding', currFolding, 'smParam', currSmParam)

    close all 
end

% mapXcorrCurvePermutation_Vel (xcorr chan vs velocity)
mapXCorrDirName = 'mapCrossCorr';      
figuresDir2 = fullfile(p.OutputDirectory, mapXCorrDirName);

for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    chanName = currChanName{i};

    mapXcorrCurvePermutation_Vel(movieData, iChan, chanName, currMaxLayer, figuresDir2, 'impute', currImpute, ...
        'movingAvgSmoothing', currMovingAvgSmoothing, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
        'WithN', currWithN, 'numPerm', currNumPerm, 'lagMax', currLagMax, 'figFlag', currFigFlag, 'fullRange', currFullRange, ...
        'parpoolNum', currParpoolNum, 'rseed', currRseed, 'h0', currH0, 'mvFrSize', currMvFrSize, 'Folding', currFolding, ...
        'topograph', currTopograph)

    close all
end

% mapXcorrCurvePermutation (xcorr ch1 v.s. ch2 or any combination of 2 channels, if nChan > 2);
% NOTE iChan1 is bigger than iChan2, e.g. iChan1 is 2, and iChan2 is 1!;
if numel(p.ChannelIndex) >= 2
    v = 1:numel(p.ChannelIndex);
    combi = nchoosek(v,2);
    [m, ~] = size(combi); % m possible combinations
    
    for i = 1:m
        iChan1 = combi(i, 2);
        iChan2 = combi(i, 1);
        chan1Name = currChanName{iChan1};
        chan2Name = currChanName{iChan2};
        
        mapXcorrCurvePermutation(movieData, iChan1, iChan2, chan1Name, chan2Name, currMaxLayer, figuresDir2, ...
            'impute', currImpute, 'parpoolNum', currParpoolNum, 'WithN', currWithN, 'numPerm', currNumPerm, ...
            'omittedWindows', currOmittedWindows, 'movingAvgSmoothing', currMovingAvgSmoothing, ...
            'subFrames', currSubFrames, 'figFlag', currFigFlag, 'lagMax', currLagMax, 'fullRange', currFullRange, ...
            'rseed', currRseed, 'topograph', currTopograph, 'Folding', currFolding)
        
        close all
    end
end
% % MLsummary_XcorrCurvesVelAcf
% summaryDirName = 'SummaryMapDDX_Xcorr';

% % Continue here, think about how to accommodate ML into the process
% MLsummary_XcorrCurvesVelAcf(ML, iChan, iChan0, currChanName, chan0Name, currMaxLayer, ...
%     mapDescriptivesDirName, mapXCorrDirName, 'lagMax0', currLagMax, ...
%     'outDirName', summaryDirName, 'timeInterval', currTimeInterval) % Here used ML to calculate interaction btw MD1 and MD2.
% close all

end