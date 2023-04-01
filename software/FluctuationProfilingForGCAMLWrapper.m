function FluctuationProfilingForGCAMLWrapper(movieListOrProcess, varargin)
% FluctuationProfilingForGCAMLWrapper is a wrapper function for phaseMasking, 
% phaseDescriptives_OneChan, phaseDescriptives_MaxMinVel_OneChan, 
% and MLsummary_FluctuationProfiling.
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
% J Noh, 7/2021. Modified for GCA package.
% As of July 2016, the first argument could also be a Process. Use
% getOwnerAndProcess to simplify compatability.
%
% As of August 2016, the standard second argument should be the parameter
% structure
%
% Qiongjing (Jenny) Zou, Sep 2018
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
    'FluctuationProfilingForGCAProcessML',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used 
% rather than the one stored in FluctuationProfilingProcessML

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
outName = 'FluctuationProfiling';       

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
    
    outPathsMD = cell(numel(movieData.channels_), 1);   % per chanID analyzed
    for i = 1:numel(p.ChannelIndex)
        outPathsMD{p.ChannelIndex(i)} = [movieData.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
                filesep ['FluctuationProfiling_', p.chanNameWithPreprocess{i}]];
        mkClrDir(outPathsMD{p.ChannelIndex(i)}); % this func does not clear subfolders!             
    end
    outFilePaths{1, j} = outPathsMD;
end

% calculations at ML level:
% logging input paths (bookkeeping) - continue
process.setInFilePaths(allinFilePaths);

% logging output paths - continue
currOutputDirectory = cell(numel(movieData.channels_), 1);
currOutputDirectory{1, 1} = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep ['FluctuationProfiling_', p.chanNameWithPreprocess{1}]];
mkClrDir(currOutputDirectory{1, 1});    
currOutputDirectory{2, 1} = [movieList.outputDirectory_ filesep 'GrangerCausalityAnalysisPackage' ...
        filesep ['FluctuationProfiling_', p.chanNameWithPreprocess{2}]];
mkClrDir(currOutputDirectory{2, 1});        
    
outFilePaths{numMDs+1} = currOutputDirectory;
 
process.setOutFilePaths(outFilePaths);

%% Algorithm

for j = 1:numMDs
    movieData = MDs{j};
    disp('============================')
    disp(['======= MD index = ', num2str(j), ' ======='])
    disp('============================')
    
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
    
    % Fluctuation Profiling for each channel
    for i = 1:numel(p.ChannelIndex)
        
        figuresDir = outFilePaths{j}{p.ChannelIndex(i)};        
        iChan = p.chanCodeWithPreprocess(i);
        chanName = p.chanNameWithPreprocess{i};    
        
        % phaseMasking to determine Prot/Ret phases

        smParamTh = p.smParamPhaseSegmentation;
        [protMask1, retMask1] = phaseMasking(movieData, smParamTh, figuresDir, ...
            'omittedWindows', currOmittedWindows, ...
            'minimumRunLength', p.minimumRunLength, ... 
            'EWMA', p.EWMAlambda, ...
            'movingAvgSmoothing', false);               % pre-fixed pars in GCA or this process 
            
        % phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask, samplingBw, figuresDir, varargin)
        
        %  Retraction  
        Mask0 = retMask1;
        chanNamePhase = [chanName, '-ret'];
        phaseDescriptives_OneChan(movieData, iChan, p.maxLayer, chanNamePhase, ...
            chanName, Mask0, p.samplingBw, figuresDir, ...
            'omittedWindows', currOmittedWindows, ...
            'impute', p.impute, ...
            'WithN', p.WithN, ...
            'figFlag', p.figFlag,  ...
            'movMedFrameSize', p.movMedFrameSize, ... 
            'EWMA', p.EWMAlambda, ...
            'movingAvgSmoothing', false)                % pre-fixed pars in GCA or this process 
        
        %  Protrusion  
        Mask0 = protMask1;
        chanNamePhase = [chanName, '-prot'];
        phaseDescriptives_OneChan(movieData, iChan, p.maxLayer, chanNamePhase, ...
            chanName, Mask0, p.samplingBw, figuresDir, ...
            'omittedWindows', currOmittedWindows, ...
            'impute', p.impute, ...
            'WithN', p.WithN, ...
            'figFlag', p.figFlag,  ...
            'movMedFrameSize', p.movMedFrameSize, ... 
            'EWMA', p.EWMAlambda, ...
            'movingAvgSmoothing', false)                % pre-fixed pars in GCA or this process 

        % Min/Maximum Velocity
        % phaseDescriptives_MaxMinVel_OneChan(MD, iChan, maxLayer, chanName, ...
        %   chanTitle, smParamTh, samplingBw, figuresDir, varargin) 
        chanNamePhase = chanName;
        phaseDescriptives_MaxMinVel_OneChan(movieData, iChan, p.maxLayer, chanNamePhase, ...
            chanName, p.smParamPhaseSegmentation, p.samplingBw, figuresDir, ...
            'minimumRunLength', p.minimumRunLength, ...
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

for i = 1:numel(p.ChannelIndex)
    % iChan = p.chanCodeWithPreprocess(i);
    chanName = p.chanNameWithPreprocess{i};
    MDDirName = ['GrangerCausalityAnalysisPackage' ...
                    filesep ['FluctuationProfiling_', p.chanNameWithPreprocess{i}]];
    
    % MLsummary_FluctuationProfiling(ML, chanName, maxLayer, analNameDesc, outDirName, varargin)                
    MLsummary_FluctuationProfiling(movieList, chanName, p.maxLayer, ...
        MDDirName, MDDirName, ...
        'lagMax0', p.samplingBw, ...
        'figFlag', p.figFlag)   
end

end