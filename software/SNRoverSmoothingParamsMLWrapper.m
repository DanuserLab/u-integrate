function SNRoverSmoothingParamsMLWrapper(movieListOrProcess, varargin)
% SNRoverSmoothingParamsMLWrapper is a wrapper for
% MD_SNRoverSmoothingParams_EWMA.m and MLsummary_SNRoverSmoothingParams.m.
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

[movieList, process] = getOwnerAndProcess(movieListOrProcess,'SNRoverSmoothingParamsProcessML',true);
p = parseProcessParams(process, paramsIn);  % If parameters are explicitly given, they should be used
                                        % rather than the one stored in XcorrAnalysisProcessML

% p.avgLengthOfProtRetCycle is required for setting up multiple pars like 
% movMedFrameSize, ARorder, lagMax, etc                                         
%if ~isfinite(p.avgLengthOfProtRetCycle)
%    error('avgLengthOfProtRetCycle (in frames, an integer) is required.')
%end
                                        
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
    inFilePaths = cell(2, nChan);   % (prot or winSampling, chanID analyzed)
    for i = p.ChannelIndex
        protSamplingProc = movieData.getProcess(protSamplingId);
        inFilePaths{1,i} = protSamplingProc.outFilePaths_;  % there is only one outFilePaths_ in ProtrusionSamplingProcess for all channels.
        
        winSamplingProc = movieData.getProcess(winSamplingId);
        inFilePaths{2,i} = winSamplingProc.outFilePaths_{1,i};
    end
    allinFilePaths{1,j} = inFilePaths;
    
    % logging output paths
    % Did not separate MD level output for channels.
    outFilePaths{j} = [movieData.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep 'SNRoverSmoothingParameters'];
    mkClrDir(outFilePaths{j}); % this func does not clear subfolders! 

end

% calculations at ML level:
% logging input paths (bookkeeping) - continue
process.setInFilePaths(allinFilePaths);

% logging output paths - continue
currOutputDirectory = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep 'SNRoverSmoothingParameters'];
    
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
    
    % input: MD, iChan, maxLayer, chanName, arOrder, figuresDir, varargin
    % ommittedWindows = [] by default in computing SNRs
    MD_SNRoverSmoothingParams_EWMA(movieData, ...
        0, 1, 'Vel', p.ARorder, figuresDir, ...
        'impute', p.impute, ...
        'WithN', p.WithN, ...
        'figFlag', p.figFlag,  ...
        'movMedFrameSize', p.movMedFrameSize, ...
        'movingAvgSmoothing', false)            % pre-fixed pars in GCA or this process 
    
    for i = 1:numel(p.ChannelIndex)
        iChan = p.chanCodeWithPreprocess(i);
        chanName = p.chanNameWithPreprocess{i};
        
        MD_SNRoverSmoothingParams_EWMA(movieData, ...
            iChan, p.maxLayer, chanName, p.ARorder, figuresDir, ...
            'impute', p.impute, ...
            'WithN', p.WithN, ...
            'figFlag', p.figFlag,  ...
            'movMedFrameSize', p.movMedFrameSize, ...
            'movingAvgSmoothing', false)            % pre-fixed pars in GCA or this process 
    end 
end

%% Summarize at the movieList level

% outDirName for both of MD and ML
outDirName = ['GrangerCausalityAnalysisPackage' ...
        filesep 'SNRoverSmoothingParameters'];

MLsummary_SNRoverSmoothingParams(movieList, 0, 'Vel', outDirName, ...
    'outDirName', outDirName)           % ML, iChan, chanName, SNRdirName, varargin

for i = 1:numel(p.ChannelIndex)
    iChan = p.chanCodeWithPreprocess(i);
    chanName = p.chanNameWithPreprocess{i};
    
    MLsummary_SNRoverSmoothingParams(movieList, iChan, chanName, outDirName, ...
        'outDirName', outDirName)
end

end

