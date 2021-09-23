function XcorrAnalysisML(movieListOrProcess, varargin)
% XcorrAnalysis wrapper function for 
% mapXcorrCurvePermutation_Vel, and MLsummary_XcorrCurvesVelAcf to be executed by
% XcorrAnalysisProcessML.
%
% INPUT
% movieListOrProcess - either a MovieList (legacy)
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
% 10/29/2018 Deleted the Xcorr ch1 v.s. ch2 (chN v.s. ChN+1, when N>=1) as
% per Jungsik's request.
%
% Qiongjing (Jenny) Zou, Aug 2018


%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('movieListOrProcess', @isProcessOrMovieList);
ip.addOptional('param',[], @isstruct);
ip.parse(movieListOrProcess, varargin{:});
p = ip.Results;
paramsIn = p.param;

%% Registration

% Get MovieList object and Process
% If movieListOrProcess is a MovieList and does not contain an
% XcorrAnalysisProcessML, create XcorrAnalysisProcessML using constructor with no
% arguments.
% If movieListOrProcess is a MovieList and does contain an XcorrAnalysisProcessML,
% then return the first instance of an XcorrAnalysisProcessML.
% If movieListOrProcess is an XcorrAnalysisProcessML, then return the Process and it's
% MovieList owner.
% Otherwise throw an error.
[movieList, process] = getOwnerAndProcess(movieListOrProcess,'XcorrAnalysisProcessML',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used
% rather than the one stored in XcorrAnalysisProcessML

% Parameters
currMaxLayer = p.maxLayer;
currImpute = p.impute;
currMovingAvgSmoothing = p.movingAvgSmoothing;
currSubFrames = p.subFrames;
currNumPerm = p.numPerm;
currTopograph = p.topograph;
currFigFlag = p.figFlag;
currWithN = p.WithN;
% currParpoolNum = p.parpoolNum; % parpoolNum is no longer a needed parameter for mapXcorrCurvePermutation_Vel
currRseed = p.rseed;
currFolding = p.Folding;
currLagMax = p.lagMax;
currFullRange = p.fullRange;
currH0 = p.h0;
currMvFrSize = p.mvFrSize;
currChanName = p.chanName; % cell array, # cell = # channels
currTimeInterval = p.timeInterval;


% load MDs from movieList:
% MDs = movieList.getMovies();
% numMDs = numel(MDs);

% if p.MovieDataIndex used, make a new MDs, and ML
MDss= movieList.getMovies();
numMDs = numel(p.MovieDataIndex);
MDs = cell(1, numMDs);
for i = 1:numMDs
    MDs{i} = MDss{p.MovieDataIndex(i)}
end
MLnew = MovieList(MDs, movieList.movieListPath_);
% error checking for p.MovieDataIndex:
activityMapId = movieList.getProcessIndex('ActivityMapDescriptionProcessML');
activityMapProc=movieList.getProcess(activityMapId);
if ~all(ismember(p.MovieDataIndex, activityMapProc.funParams_.MovieDataIndex))
    error("Please select MovieData which has finished the Activity Map Description Process.")
end

allinFilePaths = cell(1,numMDs); % prepare to log input paths (bookkeeping).

% calculations at MD level:
for j = 1:numMDs
    movieData = MDs{j};
    
    % Sanity Checks
    nChan = numel(movieData.channels_);
    if max(p.ChannelIndex) > nChan || max(p.ChannelIndex) < nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
        error('Invalid channel numbers specified! Check ChannelIndex input!!')
    end

    if numel(currChanName) ~= numel(movieData.channels_);
        error('Please provide a valid input for ''Channel Name(s)''.') % this happens b/c something wrong with user input of ChanName in GUI.
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
    inFilePaths = cell(3, nChan);
    for i = p.ChannelIndex
        protSamplingProc = movieData.getProcess(protSamplingId);
        inFilePaths{1,i} = protSamplingProc.outFilePaths_;% there is only one outFilePaths_ in ProtrusionSamplingProcess for all channels.
        
        winSamplingProc = movieData.getProcess(winSamplingId);
        inFilePaths{2,i} = winSamplingProc.outFilePaths_{1,i};
        
        LBTestId = movieList.getProcessIndex('QuiescentWindowDetectionProcessML');
        if ~isempty(LBTestId)
            LBTestProc = movieList.getProcess(LBTestId);
            if LBTestProc.success_ == 1
                inFilePaths{3,i} = LBTestProc.outFilePaths_{j};% there is only one outFilePaths_ in LBTestProcess for all channels.
            else
                inFilePaths{3,i} = [];
            end
        else
            inFilePaths{3,i} = [];
        end
    end
    allinFilePaths{1,j} = inFilePaths;
    
    % logging output paths
    % Did not separate MD level output for channels.
    outFilePaths{j} = [movieData.outputDirectory_ filesep 'XcorrFluctuationPackage' filesep 'XcorrAnalysis'];
    mkClrDir(outFilePaths{j}); % this func does not clear subfolders! 
    
    if p.omittedWindows == false
        currOmittedWindows = [];
    else
        % If LBTestProcess (QuiescentWindowDetectionProcessML) is done, exclude quiescent windows (calculate currMovingAvgSmoothing).
        if ~isempty(LBTestId)
            LBTestProc = movieList.getProcess(LBTestId);
            if LBTestProc.success_ == 1
                LBTestOutDir = LBTestProc.outFilePaths_{j};
                indPath = fullfile(LBTestOutDir, 'indActive_windowIndex.mat');
                raw = load(indPath);
                
                omitWin = find(raw.indActive == 0);
                disp('== omittedWindows: ==')
                disp(omitWin)
                currOmittedWindows = omitWin;
            else
                error("Quiescent Window Detection needs to be successfully finished to exclude quiescent windows in the Cross Correlation Analysis.")
            end
        else
            error("Quiescent Window Detection needs to be successfully finished to exclude quiescent windows in the Cross Correlation Analysis.")
        end
    end

    %% Algorithm
    
    % mapXcorrCurvePermutation_Vel (xcorr chan vs velocity)
    figuresDir2 = outFilePaths{j};
    
    for i = 1:numel(p.ChannelIndex)
        iChan = p.ChannelIndex(i);
        chanName = currChanName{i};
        
        mapXcorrCurvePermutation_Vel(movieData, iChan, chanName, currMaxLayer, figuresDir2, 'impute', currImpute, ...
            'movingAvgSmoothing', currMovingAvgSmoothing, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
            'WithN', currWithN, 'numPerm', currNumPerm, 'lagMax', currLagMax, 'figFlag', currFigFlag, 'fullRange', currFullRange, ...
            'rseed', currRseed, 'h0', currH0, 'mvFrSize', currMvFrSize, 'Folding', currFolding, ...
            'topograph', currTopograph)
        
        % close all % this will close GUI window, so comment out.
    end
    
    %%% Commented below to delete the Xcorr ch1 v.s. ch2 (chN v.s. ChN+1, when N>=1)
    % % mapXcorrCurvePermutation (xcorr ch1 v.s. ch2 or any combination of 2 channels, if nChan >= 2);
    % % NOTE iChan1 is bigger than iChan2, e.g. iChan1 is 2, and iChan2 is 1!;
    
    % if numel(p.ChannelIndex) >= 2
    %     v = 1:numel(p.ChannelIndex);
    %     combi = nchoosek(v,2);
    %     [m, ~] = size(combi); % m possible combinations
        
    %     for i = 1:m
    %         iChan1 = combi(i, 2);
    %         iChan2 = combi(i, 1);
    %         chan1Name = currChanName{iChan1};
    %         chan2Name = currChanName{iChan2};
            
    %         mapXcorrCurvePermutation(movieData, iChan1, iChan2, chan1Name, chan2Name, currMaxLayer, figuresDir2, ...
    %             'impute', currImpute, 'parpoolNum', currParpoolNum, 'WithN', currWithN, 'numPerm', currNumPerm, ...
    %             'omittedWindows', currOmittedWindows, 'movingAvgSmoothing', currMovingAvgSmoothing, ...
    %             'subFrames', currSubFrames, 'figFlag', currFigFlag, 'lagMax', currLagMax, 'fullRange', currFullRange, ...
    %             'rseed', currRseed, 'topograph', currTopograph, 'Folding', currFolding)
            
    %         % close all % this will close GUI window, so comment out.
    %     end
    % end
    
end

% calculations at ML level:

% logging input paths (bookkeeping) - continue
process.setInFilePaths(allinFilePaths);

% logging output paths - continue
for i = 1:nChan
    currOutputDirectory{i} = [p.OutputDirectory  filesep 'Summary_Xcorr_ch' num2str(i) 'ch0'];
end

for k = numMDs+1:numMDs+nChan
    outFilePaths{k} = currOutputDirectory{k-numMDs};
    mkClrDir(outFilePaths{k});
end
process.setOutFilePaths(outFilePaths);

%% Algorithm - continue

% MLsummary_XcorrCurvesVelAcf, run all combinations between 2 channels.
for i = 1:nChan
    endout=regexp(currOutputDirectory{i},filesep,'split');
    if isequal(endout{end-1}, 'XcorrFluctuationPackage')
        summaryDirName = [endout{end-1} filesep endout{end}];
    else
        summaryDirName = endout{end};
    end
    iChan = i;
    chanName = currChanName{iChan};
    iChan0 = 0;
    chan0Name = 'Vel';
    
    analNameAcf = ['XcorrFluctuationPackage' filesep 'mapDescriptives'];
    analNameXcf = ['XcorrFluctuationPackage' filesep 'XcorrAnalysis'];
    MLsummary_XcorrCurvesVelAcf(MLnew, iChan, iChan0, chanName, chan0Name, currMaxLayer, ...
        analNameAcf, analNameXcf, 'lagMax0', currLagMax, ...
        'outDirName', summaryDirName, 'timeInterval', currTimeInterval, 'figFlag', currFigFlag) % Here used ML to calculate interaction btw MD1 and MD2.
    % close all  % this will close GUI window, so comment out.
end

%%% Commented out below and changed accordingly above to delete the Xcorr ch1 v.s. ch2 (chN v.s. ChN+1, when N>=1)
% % logging output paths - continue
% v = 1: nChan+1; % +1 to include ch0
% combi = nchoosek(v,2);
% [m, ~] = size(combi); % m possible combinations
% for i = 1:m
%     currOutputDirectory{i} = [p.OutputDirectory  filesep 'Summary_Xcorr_ch' num2str(combi(i,2)-1) 'ch' num2str(combi(i,1)-1)];
% end

% for k = numMDs+1:numMDs+m
%     outFilePaths{k} = currOutputDirectory{k-numMDs};
%     mkClrDir(outFilePaths{k});
% end
% process.setOutFilePaths(outFilePaths);

% %% Algorithm - continue

% % MLsummary_XcorrCurvesVelAcf, run all combinations between 2 channels.
% for i = 1:m
%     endout=regexp(currOutputDirectory{i},filesep,'split');
%     if isequal(endout{end-1}, 'XcorrFluctuationPackage')
%         summaryDirName = [endout{end-1} filesep endout{end}];
%     else
%         summaryDirName = endout{end};
%     end
%     iChan = combi(i,2)-1;
%     chanName = currChanName{iChan};
%     iChan0 = combi(i,1)-1;
%     if iChan0 == 0
%         chan0Name = 'Vel';
%     else
%         chan0Name = currChanName{iChan0};
%     end
    
%     analNameAcf = ['XcorrFluctuationPackage' filesep 'Xcorr_Analysis' filesep mapDescriptivesDirName];
%     analNameXcf = ['XcorrFluctuationPackage' filesep 'Xcorr_Analysis' filesep mapXCorrDirName];
%     MLsummary_XcorrCurvesVelAcf(movieList, iChan, iChan0, chanName, chan0Name, currMaxLayer, ...
%         analNameAcf, analNameXcf, 'lagMax0', currLagMax, ...
%         'outDirName', summaryDirName, 'timeInterval', currTimeInterval, 'figFlag', currFigFlag) % Here used ML to calculate interaction btw MD1 and MD2.
%     % close all  % this will close GUI window, so comment out.
% end

end