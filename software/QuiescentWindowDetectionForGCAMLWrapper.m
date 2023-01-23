function QuiescentWindowDetectionForGCAMLWrapper(movieListOrProcess, varargin)
% QuiescentWindowDetectionForGCAMLWrapper is a wrapper function for 
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
% J Noh, 7/2021. Modify it for GCA package.
% Qiongjing (Jenny) Zou, Aug 2018
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
                            'QuiescentWindowDetectionForGCAProcessML',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used 
                % rather than the one stored in QuiescentWindowDetectionProcessML

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
outName = 'QuiescentWindowDetection';

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
    inFilePaths = cell(2, numel(p.ChannelIndex));   % (prot or winSampling, chanID analyzed)
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

for j = 1:numMDs
    movieData = MDs{j};
    disp('============================')
    disp(['======= MD index = ', num2str(j), ' ======='])
    disp('============================')
    
    figuresDir = outFilePaths{j};
    
    % mapDescriptives_Vel_LB(MD, figuresDir, varargin)
    mapDescriptives_Vel_LB(movieData, figuresDir, ...
        'impute', p.impute, ...  
        'figFlag', p.figFlag,  ...
        'movmeanNum', p.movmeanNum, ...         % # of windows for moving averages of P-values for smooth outcome
        'topograph', p.topograph, ...
        'EWMA', p.EWMAlambda, ...
        'movingAvgSmoothing', false)            % pre-fixed pars in GCA or this process 
end

%% Summarize at the movieList level

% outDirName for both of MD and ML
outDirName = [filesep 'GrangerCausalityAnalysisPackage' ...
        filesep outName];
    
% MLsummary_quiescentWindow(ML, LBdirName, varargin)        
MLsummary_quiescentWindow(movieList, outDirName, 'outDirName', outDirName)

end