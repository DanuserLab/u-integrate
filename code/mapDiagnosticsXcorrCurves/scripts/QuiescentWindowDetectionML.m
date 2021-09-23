function QuiescentWindowDetectionML(movieListOrProcess, varargin)
% QuiescentWindowDetectionML wrapper function for mapDescriptives_Vel_LB to be executed by
% QuiescentWindowDetectionProcessML.
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
% QuiescentWindowDetectionML.m was created to replace LBTestML.m on 10/24/2018.
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
% QuiescentWindowDetectionProcessML, create QuiescentWindowDetectionProcessML using constructor with no
% arguments.
% If movieListOrProcess is a MovieList and does contain an QuiescentWindowDetectionProcessML,
% then return the first instance of an QuiescentWindowDetectionProcessML.
% If movieListOrProcess is an QuiescentWindowDetectionProcessML, then return the Process and it's
% MovieList owner.
% Otherwise throw an error.
[movieList, process] = getOwnerAndProcess(movieListOrProcess,'QuiescentWindowDetectionProcessML',true);
p = parseProcessParams(process, paramsIn); % If parameters are explicitly given, they should be used 
% rather than the one stored in QuiescentWindowDetectionProcessML

% Parameters
currImpute = p.impute;
currMovingAvgSmoothing = p.movingAvgSmoothing;
currFigFlag = p.figFlag;
currOmittedWindows = p.omittedWindows;
currSubFrames = p.subFrames;
currTopograph = p.topograph;
currFolding = p.Folding;
currDerivative = p.derivative;
currSmParam = p.smParam;
currMovmeanNum = p.movmeanNum;
currRatioThreshold = p.ratioThreshold;

% % Sanity Checks
% nChan = numel(movieData.channels_);
% if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
%     error('Invalid channel numbers specified! Check ChannelIndex input!!')
% end

% load MDs from movieList:
MDs = movieList.getMovies();
numMDs = numel(MDs);

for i = 1:numMDs
	movieData = MDs{i};

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

	if isempty(movieData.pixelSize_); error('movieData.pixelSize_ is required.'); end
	if isempty(movieData.timeInterval_); error('movieData.timeInterval_ is required.'); end

	% logging input paths (bookkeeping)
	protSamplingProc = movieData.getProcess(protSamplingId);
	inFilePaths{i} = protSamplingProc.outFilePaths_;% there is only one outFilePaths_ in ProtrusionSamplingProcess for all channels.
	
	% logging output paths
	outFilePaths{i} = [movieData.outputDirectory_ filesep 'XcorrFluctuationPackage' filesep 'QuiescentWindowDetection']; 
	mkClrDir(outFilePaths{i});
	

	%% Algorithm
	figuresDir = outFilePaths{i};
	mapDescriptives_Vel_LB(movieData, figuresDir, 'impute', currImpute,'movingAvgSmoothing', currMovingAvgSmoothing, ...
	        'figFlag', currFigFlag, 'omittedWindows', currOmittedWindows, 'subFrames', currSubFrames, ...
	        'topograph', currTopograph, 'derivative', currDerivative, 'Folding', currFolding, 'smParam', currSmParam, ...
	        'movmeanNum', currMovmeanNum)

end

% logging input paths (bookkeeping) - continue
process.setInFilePaths(inFilePaths);

% logging output paths - continue
currOutputDirectory = p.OutputDirectory;
outFilePaths{numMDs+1} = [p.OutputDirectory filesep 'Summary_QuiescentWindow'];
mkClrDir(outFilePaths{numMDs+1});
process.setOutFilePaths(outFilePaths);

%% Algorithm - continue

% MLsummary_quiescentWindow
LBdirName = ['XcorrFluctuationPackage' filesep 'QuiescentWindowDetection'];
endout=regexp(outFilePaths{numMDs+1},filesep,'split');
if isequal(endout{end-1}, 'XcorrFluctuationPackage')
    summaryDirName = [endout{end-1} filesep endout{end}];
else
    summaryDirName = endout{end};
end
MLsummary_quiescentWindow(movieList, LBdirName, 'ratioThreshold', currRatioThreshold, 'outDirName', summaryDirName)

% Not nessary to put differnet channel results in different folder as below, b/c all channel will 
% have the same QuiescentWindowDetection results.

% dName = 'QuiescentWindowDetection_for_channel_';
% nChan = numel(movieData.channels_);
% outFilePaths = cell(1, nChan);
% for iChan = p.ChannelIndex
% 	currDir = [p.OutputDirectory filesep dName num2str(iChan)];
%     outFilePaths{1,iChan} = currDir;
%     process.setOutFilePaths(currDir, iChan);
%     mkClrDir(outFilePaths{1,iChan});
%     % Output will be saved to currDir/indActive_windowIndex.mat
% 	% outputFile = fullfile(currDir, 'indActive_windowIndex.mat'); % output

% %% Algorithm
% figuresDir = currDir; % output dir for this process

% mapDescriptives_Vel_LB(movieData, figuresDir, 'impute',1,'movingAvgSmoothing',1, ...
%         'figFlag','off', 'omittedWindows', [], 'subFrames', [], 'topograph', 'on')
% end

end