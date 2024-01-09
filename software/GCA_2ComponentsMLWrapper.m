function GCA_2ComponentsMLWrapper(movieListOrProcess, varargin)
% GCA_3ComponentsMLWrapper is a wrapper for MD_iGC_SPAR3ch.m and
% MLsummary_iGC_SPAR3ch.m.
%
%
% Jungsik Noh, 7/2021
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
%
% This file is part of GrangerCausalityAnalysisPackage.
% 
% GrangerCausalityAnalysisPackage is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GrangerCausalityAnalysisPackage is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GrangerCausalityAnalysisPackage.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
 
%% Input
ip = inputParser;
ip.addRequired('movieListOrProcess', @isProcessOrMovieList);
ip.addOptional('param',[], @isstruct);
ip.parse(movieListOrProcess, varargin{:});
p = ip.Results;
paramsIn = p.param;

%% Registration

[movieList, process] = getOwnerAndProcess(movieListOrProcess,'GCA_2ComponentsProcessML',true);
p = parseProcessParams(process, paramsIn);  % If parameters are explicitly given, they should be used
                                        % rather than the one stored in XcorrAnalysisProcessML
                                        
LBTestId = movieList.getProcessIndex('QuiescentWindowDetectionForGCAProcessML');
LBTestProc = [];
if ~isempty(LBTestId)
    LBTestProc = movieList.getProcess(LBTestId);
end  

% (obsolete) if p.MovieDataIndex used, make a new MDs, and ML
MDss= movieList.getMovies();
numMDs = numel(p.MovieDataIndex);
MDs = cell(1, numMDs);
for i = 1:numMDs
    MDs{i} = MDss{p.MovieDataIndex(i)};
end
%MLnew = MovieList(MDs, movieList.movieListPath_);

allinFilePaths = cell(1,numMDs); % prepare to log input paths (bookkeeping).
outFilePaths = cell(1, numMDs + 1);
outName = ['GCA_2Components_', p.chanNameWithPreprocess{1}, '_', p.chanNameWithPreprocess{2}];

% checking at MD level:
for j = 1:numMDs
    movieData = MDs{j};
    
    % Sanity Checks
    disp(['== Checking movieData: ', num2str(j)])
    
    nChan = numel(movieData.channels_);
    chans = 1:nChan;
    if ~isempty(setdiff(p.ChannelIndex, chans))
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
    inFilePaths = cell(3, numel(p.ChannelIndex));   % (prot or winSampling, chanID analyzed)
    for i = p.ChannelIndex
        protSamplingProc = movieData.getProcess(protSamplingId);
        inFilePaths{1,i} = protSamplingProc.outFilePaths_;  % there is only one outFilePaths_ in ProtrusionSamplingProcess for all channels.
        
        winSamplingProc = movieData.getProcess(winSamplingId);
        inFilePaths{2,i} = winSamplingProc.outFilePaths_{1,i};
        
        if ~isempty(LBTestProc)
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
    outFilePaths{j} = [movieData.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep outName];
    mkClrDir(outFilePaths{j}); % this func does not clear subfolders! 

end

% calculations at ML level:
% logging input paths (bookkeeping) - continue
process.setInFilePaths(allinFilePaths);

% logging output paths - continue
currOutputDirectory = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep outName];
    
outFilePaths{numMDs+1} = currOutputDirectory;
mkClrDir(outFilePaths{numMDs+1});
process.setOutFilePaths(outFilePaths);

%% Algorithm

chIds = [p.chanCodeWithPreprocess, 0];
chNames = [p.chanNameWithPreprocess, {'Vel'}];
relationMat =  [2, 1; 
                1, 2;
                1, 3;
                3, 1;
                2, 3;
                3, 2];

for j = 1:numMDs
    movieData = MDs{j};
    disp('============================')
    disp(['======= MD index = ', num2str(j), ' ======='])
    disp('============================')
    
    figuresDir = outFilePaths{j};
    
    if p.excludeQuiescentWindows == false
        currOmittedWindows = [];
    else 
        if ~isempty(LBTestProc)
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
        
    for r = 1:size(relationMat, 1)
        id1 = relationMat(r, 1);
        id2 = relationMat(r, 2);
        %id3 = relationMat(r, 3);
        iChan1 = chIds(id1); iChan2 = chIds(id2); %iChan3 = chIds(id3);   
        chan1Name = chNames{id1}; chan2Name = chNames{id2}; %chan3Name = chNames{id3};
        
        timeWindowLagMax = [p.SPARLagMax, 1];   % [maximum time lag for AR, window lag]

        % MD_iGC_SPAR2ch(MD, iChan1, iChan2,  chan1Name, chan2Name, layerMax, ...
        % figuresDir, twlagMax, twlagMaxReg, varargin)
        % GC from iChan1 to iChan2 given iChan3
        MD_iGC_SPAR2ch(movieData, iChan1, iChan2, chan1Name, chan2Name, ...
            p.maxLayer, figuresDir, timeWindowLagMax, timeWindowLagMax, ...
            'omittedWindows', currOmittedWindows, ...
            'impute', p.impute, ...
            'WithN', p.WithN, ...
            'figFlag', p.figFlag,  ...
            'movMedFrameSize', p.movMedFrameSize, ...
            'EWMA', p.EWMAlambda, ...
            'movingAvgSmoothing', false)                % pre-fixed pars in GCA or this process 
        
    end 
end

%% Summarize at the movieList level

% outDirName for both of MD and ML
outDirName = [filesep 'GrangerCausalityAnalysisPackage' ...
        filesep outName];

for r = 1:size(relationMat, 1)
    id1 = relationMat(r, 1);
    id2 = relationMat(r, 2);
    %id3 = relationMat(r, 3);
    chan1Name = chNames{id1}; chan2Name = chNames{id2}; %chan3Name = chNames{id3};

    MLsummary_iGC_SPAR2ch(movieList, chan1Name, chan2Name, ...
        p.maxLayer, outDirName,'outDirName', outDirName)  
end

%% Draw Granger-causality network diagrams

disp('==== drawing network diagram')

chanDetailedNames = chNames;
chanNamesforNet = [p.ChannelName, 'Edge Velocity'];
GCsummaryDir = outFilePaths{numMDs+1};

GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, p.maxLayer)

end