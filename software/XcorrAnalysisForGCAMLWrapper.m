function XcorrAnalysisForGCAMLWrapper(movieListOrProcess, varargin)
% XcorrAnalysis wrapper function is a wrapper function for
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
% Updates:
% J Noh, 7/2021. Modified for GCA package.
% Qiongjing (Jenny) Zou, Aug 2018
%
% Copyright (C) 2022, Danuser Lab - UTSouthwestern 
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


%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('movieListOrProcess', @isProcessOrMovieList);
ip.addOptional('param',[], @isstruct);
ip.parse(movieListOrProcess, varargin{:});
p = ip.Results;
paramsIn = p.param;

%% Registration

[movieList, process] = getOwnerAndProcess(movieListOrProcess, ...
                                        'XcorrAnalysisForGCAProcessML',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used
                % rather than the one stored in XcorrAnalysisProcessML

LBTestId = movieList.getProcessIndex('QuiescentWindowDetectionForGCAProcessML');
LBTestProc = [];
if ~isempty(LBTestId)
    LBTestProc = movieList.getProcess(LBTestId);
end  

% error checking  
activityMapId = movieList.getProcessIndex('ActivityMapDescriptionForGCAProcessML');
activityMapProc = movieList.getProcess(activityMapId);
if ~all(ismember(p.MovieDataIndex, activityMapProc.funParams_.MovieDataIndex))
    error("Please select MovieData which has finished the Activity Map Description Process.")
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
outName = 'XcorrAnalysis';       

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
    inFilePaths = cell(3, numel(p.ChannelIndex));   % (prot, winSampling, LBTestProc, chanID analyzed)
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
currOutputDirectory{1} = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep ['XcorrAnalysis_', p.chanNameWithPreprocess{1}, '_Vel']];
currOutputDirectory{2} = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep ['XcorrAnalysis_', p.chanNameWithPreprocess{2}, '_Vel']];
currOutputDirectory{3} = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep ['XcorrAnalysis_', p.chanNameWithPreprocess{2}, '_', p.chanNameWithPreprocess{1}]];    
    
outFilePaths{numMDs+1} = currOutputDirectory{1};
outFilePaths{numMDs+2} = currOutputDirectory{2};
outFilePaths{numMDs+3} = currOutputDirectory{3};

mkClrDir(outFilePaths{numMDs+1});
process.setOutFilePaths(outFilePaths);

%% Algorithm

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
    
    % Xcorr with Vel
    for i = 1:numel(p.ChannelIndex)
        iChan = p.chanCodeWithPreprocess(i);
        chanName = p.chanNameWithPreprocess{i};
        
        % mapXcorrCurvePermutation_Vel(MD, iChan1, chan1Name, layerMax, figuresDir, varargin) 
        mapXcorrCurvePermutation_Vel(movieData, iChan, chanName, p.maxLayer, figuresDir, ...
            'omittedWindows', currOmittedWindows, ...
            'impute', p.impute, ...
            'WithN', p.WithN, ...
            'figFlag', p.figFlag,  ...
            'movMedFrameSize', p.movMedFrameSize, ...
            'topograph', p.topograph, ...
            'numPerm', p.numPerm, ...
            'lagMax', p.lagMax, ...
            'EWMA', p.EWMAlambda, ...
            'movingAvgSmoothing', false)                % pre-fixed pars in GCA or this process 
    end 
    
    % Xcorr between ch1/ch2
    % mapXcorrCurvePermutation(MD, iChan1, iChan2, chan1Name, chan2Name, layerMax, figuresDir, varargin) 
    iChan1 = p.chanCodeWithPreprocess(2);
    iChan2 = p.chanCodeWithPreprocess(1);
    chan1Name = p.chanNameWithPreprocess{2};
    chan2Name = p.chanNameWithPreprocess{1};
    
    mapXcorrCurvePermutation(movieData, iChan1, iChan2, chan1Name, chan2Name, p.maxLayer, figuresDir, ...
        'omittedWindows', currOmittedWindows, ...
        'impute', p.impute, ...
        'WithN', p.WithN, ...
        'figFlag', p.figFlag,  ...
        'movMedFrameSize', p.movMedFrameSize, ...
        'topograph', p.topograph, ...
        'numPerm', p.numPerm, ...
        'lagMax', p.lagMax, ...
        'EWMA', p.EWMAlambda, ...
        'movingAvgSmoothing', false)                % pre-fixed pars in GCA or this process 
    
end

%% Summarize at the movieList level

% outDirName for both of MD and ML
tmp = split(activityMapProc.outFilePaths_{1}, filesep);
activityMapProc_MDoutName = tmp{end};
analNameAcf = ['GrangerCausalityAnalysisPackage' filesep activityMapProc_MDoutName];
analNameXcf = ['GrangerCausalityAnalysisPackage' filesep 'XcorrAnalysis'];

% Xcorr with Vel
for i = 1:numel(p.ChannelIndex)
    iChan = p.chanCodeWithPreprocess(i);
    chanName = p.chanNameWithPreprocess{i};
    summaryDirName = ['GrangerCausalityAnalysisPackage' ...
                    filesep ['XcorrAnalysis_', p.chanNameWithPreprocess{i}, '_Vel']];
    
    % MLsummary_XcorrCurvesVelAcf(ML, iChan1, iChan2, chan1Name, chan2Name, ...
    % maxLayer, analNameAcf, analNameXcf, varargin)        
    MLsummary_XcorrCurvesVelAcf(movieList, iChan, 0, chanName, 'Vel', ...
        p.maxLayer, analNameAcf, analNameXcf, ...
        'lagMax0', p.lagMax, ...
        'outDirName', summaryDirName, ...
        'figFlag', p.figFlag)
        
end

% Xcorr between ch1/ch2
iChan1 = p.chanCodeWithPreprocess(2);
iChan2 = p.chanCodeWithPreprocess(1);
chan1Name = p.chanNameWithPreprocess{2};
chan2Name = p.chanNameWithPreprocess{1};
summaryDirName = ['GrangerCausalityAnalysisPackage' ...
        filesep ['XcorrAnalysis_', p.chanNameWithPreprocess{2}, '_', p.chanNameWithPreprocess{1}]];

MLsummary_XcorrCurvesVelAcf(movieList, iChan1, iChan2, chan1Name, chan2Name, ...
    p.maxLayer, analNameAcf, analNameXcf, ...
    'lagMax0', p.lagMax, ...
    'outDirName', summaryDirName, ...
    'figFlag', p.figFlag)                
 
end