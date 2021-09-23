function ActivityMapDescriptionML(movieListOrProcess, varargin)
% ActivityMapDescriptionML wrapper function for mapDescriptives_OneChan 
% to be executed by ActivityMapDescriptionProcessML.
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
% Qiongjing (Jenny) Zou, Oct 2018


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
% ActivityMapDescriptionProcessML, create ActivityMapDescriptionProcessML using constructor with no
% arguments.
% If movieListOrProcess is a MovieList and does contain an ActivityMapDescriptionProcessML,
% then return the first instance of an ActivityMapDescriptionProcessML.
% If movieListOrProcess is an ActivityMapDescriptionProcessML, then return the Process and it's
% MovieList owner.
% Otherwise throw an error.
[movieList, process] = getOwnerAndProcess(movieListOrProcess,'ActivityMapDescriptionProcessML',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used
% rather than the one stored in ActivityMapDescriptionProcessML

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
% currParpoolNum = p.parpoolNum; % parpoolNum is no longer a needed parameter for mapDescriptives_OneChan
currRseed = p.rseed;
currFolding = p.Folding;
currSmParam = p.smParam;
currChanName = p.chanName; % cell array, # cell = # channels


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

allinFilePaths = cell(1,numMDs); % prepare to log input paths (bookkeeping).

% calculations at MD level:
for j = 1:numMDs
    movieData = MDs{j};
    
    % Sanity Checks
    nChan = numel(movieData.channels_);
    if max(p.ChannelIndex) > nChan || max(p.ChannelIndex) < nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
        error('Invalid channel numbers specified! Check ChannelIndex input!!')
    end

    if numel(currChanName) ~= numel(movieData.channels_)
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
    outFilePaths{j} = [movieData.outputDirectory_ filesep 'XcorrFluctuationPackage' filesep 'mapDescriptives'];
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
                error("Quiescent Window Detection needs to be successfully finished to exclude quiescent windows in the Activity Map Description.")
            end
        else
            error("Quiescent Window Detection needs to be successfully finished to exclude quiescent windows in the Activity Map Description.")
        end
    end

    %% Algorithm
    
    % mapDescriptives_OneChan:
    figuresDir1 = outFilePaths{j};
    
    iChan0 = 0;
    maxLayer = 1;
    chan0Name = 'Vel';
    chan0Title = 'Velocity (nm/sec)';
    mapDescriptives_OneChan(movieData, iChan0, maxLayer, chan0Name, chan0Title, figuresDir1, 'adf', currAdf,'impute', currImpute, ...
        'movingAvgSmoothing', currMovingAvgSmoothing, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
        'numPerm', currNumPerm, 'topograph', currTopograph, 'figFlag', currFigFlag, 'WithN', currWithN, ...
        'rseed', currRseed, 'Folding', currFolding, 'smParam', currSmParam)
    
    
    for i = 1:numel(p.ChannelIndex)
        iChan = p.ChannelIndex(i);
        chanName = currChanName{i};
        chanTitle = chanName;
        
        mapDescriptives_OneChan(movieData, iChan, currMaxLayer, chanName, chanTitle, figuresDir1, 'adf', currAdf,'impute', currImpute, ...
            'movingAvgSmoothing', currMovingAvgSmoothing, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
            'numPerm', currNumPerm, 'topograph', currTopograph, 'figFlag', currFigFlag, 'WithN', currWithN, ...
            'rseed', currRseed, 'Folding', currFolding, 'smParam', currSmParam)
        
        % close all
    end
end

% logging input paths (bookkeeping) - continue
process.setInFilePaths(allinFilePaths);

% logging output paths - continue
outFilePaths{numMDs+1} = [p.OutputDirectory filesep 'Summary_ADFtest'];
mkClrDir(outFilePaths{numMDs+1});
process.setOutFilePaths(outFilePaths);

%% Algorithm - continue

%% MLsummary_ADFtest
mapDescDirName = ['XcorrFluctuationPackage' filesep 'mapDescriptives'];
for i = 1:nChan
    endout=regexp(outFilePaths{numMDs+1},filesep,'split');
    if isequal(endout{end-1}, 'XcorrFluctuationPackage')
        summaryDirName = [endout{end-1} filesep endout{end}];
    else
        summaryDirName = endout{end};
    end
    iChan = i;
    chanName = currChanName{iChan};

    MLsummary_ADFtest(MLnew, mapDescDirName, iChan, chanName, 'outDirName', summaryDirName)
end
end