function ActivityMapDescriptionForGCAMLWrapper(movieListOrProcess, varargin)
% ActivityMapDescriptionForGCAMLWrapper is a wrapper function for 
% mapDescriptives_OneChan.m and MLsummary_ADFtest.m.  
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
% Updates:
% J Noh, 7/2021. Modified from ActivityMapDescriptionML.m.
% Qiongjing (Jenny) Zou, Oct 2018
%
% Copyright (C) 2023, Danuser Lab - UTSouthwestern 
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
    'ActivityMapDescriptionForGCAProcessML',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used
                    % rather than the one stored in ActivityMapDescriptionProcessML
                    
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
outName = 'activityMapDescription';       

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
currOutputDirectory = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep 'Summary_ADFtest'];
    
outFilePaths{numMDs+1} = currOutputDirectory;
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
    
    % Velocity map description 
    chan0Title = 'Velocity (nm/sec)';
    
    % mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, varargin)
    mapDescriptives_OneChan(movieData, 0, 1, 'Vel', chan0Title, figuresDir, ...
        'omittedWindows', currOmittedWindows, ...
        'impute', p.impute, ...
        'WithN', p.WithN, ...
        'figFlag', p.figFlag,  ...
        'movMedFrameSize', p.movMedFrameSize, ...
        'topograph', p.topograph, ...
        'numPerm', p.numPerm, ...
        'smParam', p.smParam, ...
        'adf', p.adf, ...    
        'EWMA', p.EWMAlambda, ...
        'movingAvgSmoothing', false)                % pre-fixed pars in GCA or this process 
    
    % ch1, ch2 map description 
    for i = 1:numel(p.ChannelIndex)
        iChan = p.chanCodeWithPreprocess(i);
        chanName = p.chanNameWithPreprocess{i};
        
        mapDescriptives_OneChan(movieData, iChan, p.maxLayer, chanName, chanName, figuresDir, ...
            'omittedWindows', currOmittedWindows, ...
            'impute', p.impute, ...
            'WithN', p.WithN, ...
            'figFlag', p.figFlag,  ...
            'movMedFrameSize', p.movMedFrameSize, ...
            'topograph', p.topograph, ...
            'numPerm', p.numPerm, ...
            'smParam', p.smParam, ...
            'adf', p.adf, ...    
            'EWMA', p.EWMAlambda, ...
            'movingAvgSmoothing', false)                % pre-fixed pars in GCA or this process 
    end 
end

%% Summarize at the movieList level

% outDirName for both of MD and ML
MDoutDirName = [filesep 'GrangerCausalityAnalysisPackage' ...
        filesep 'activityMapDescription']
MLoutDirName = [filesep 'GrangerCausalityAnalysisPackage' ...
        filesep 'Summary_ADFtest'];

% MLsummary_ADFtest(ML, mapDescDirName, iChan, chanName, varargin)
MLsummary_ADFtest(movieList, MDoutDirName, 0, 'Velocity', 'outDirName', MLoutDirName)
        
for i = 1:numel(p.ChannelIndex)
    iChan = p.chanCodeWithPreprocess(i);
    chanName = p.chanNameWithPreprocess{i};
    
    MLsummary_ADFtest(movieList, MDoutDirName, iChan, chanName, 'outDirName', MLoutDirName)
end

end