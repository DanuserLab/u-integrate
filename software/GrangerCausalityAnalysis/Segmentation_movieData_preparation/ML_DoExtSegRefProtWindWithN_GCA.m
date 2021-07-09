function ML_DoExtSegRefProtWindWithN_GCA(ML, MDindex)
% MLdream1st_DoMSASegProtWindWithN IMPLEMENT fully automatic
% seg/windowing!!!
% 2017/06/22, Jungsik Noh


    disp(ML.movieListPath_) 
    disp(['MDindex:', num2str(MDindex)])
    disp(ML.movieDataFile_{MDindex})
    
    load(ML.movieDataFile_{MDindex})    % movieData.mat -> MD

% for i = 1:numel(ML.movieDataFile_);ML_DoExtSegRefProtWindWithN_GCA(ML,i); end

%% Input Parameters

iChan = 1
% and a few others below
refClRadius = 5

windparasize = 6  %6
windperpsize = 6   %4
windstcontour = 2
windsamchindex = [1;2];  % [1;2;3;4];

% Run MSA and ExternalSegmentationProcess
% Run MaskRefinementProcess


%% Run ExternalSegmentationProcess
%MSA_Seg(MD, iChan, 'tightness', msa_tightness, 'imagesOut', 1)

extSeg = ExternalSegmentationProcess(MD);
MD.addProcess(extSeg);

extSegParams = extSeg.funParams_;
extSegParams.ChannelIndex = iChan;
extSegParams.InputData = ...
{fullfile(MD.getChannelPaths{1}, '../MSAmSeg_ActinPlusmDia1cr', 'tremblingCorrected_masks131R5')};
%{fullfile(MD.outputDirectory_, 'MultiScaleAutoSeg', 'MSASeg_refined_masks_channel_1')};

% Resave the parameters
parseProcessParams(MD.processes_{end}, extSegParams);

% Run the process
MD.processes_{end}.run();
MD.sanityCheck()


%% Run MaskRefinementProcess
refineProc = MaskRefinementProcess(MD);
MD.addProcess(refineProc);

refineParams = refineProc.funParams_;

%
refineParams.ChannelIndex = iChan;
refineParams.SegProcessIndex = MD.getProcessIndex('ExternalSegmentationProcess');
refineParams.MinimumSize = 1000;
refineParams.ClosureRadius = refClRadius;
refineParams.ObjectNumber = 1;
refineParams.FillHoles = 1;

% Resave the parameterse
parseProcessParams(MD.processes_{end}, refineParams);

% Run the process
MD.processes_{end}.run();
MD.sanityCheck()


%% Run the protrusion vectors
% Create a segmentation package
prot = ProtrusionProcess(MD);
MD.addProcess(prot);

% Get the Protrusion Parameters so you can modify them
%protParams = MD.processes_{3}.funParams_;
protParams = prot.funParams_;

% Set the Protrusion Parameters
protParams.SegProcessIndex = MD.getProcessIndex('MaskRefinementProcess'); %
protParams.ChannelIndex = iChan;



% Resave the parameters
parseProcessParams(MD.processes_{end}, protParams);

% Run the process
MD.processes_{end}.run();
MD.sanityCheck()

%% Run the Windowing

wind = WindowingProcess(MD);
MD.addProcess(wind);
% % Associate the threshold process to the package
%MD.packages_{2}.setProcess(2, wind);

% Get the Windowing Parameters so you can modify them
windParams = wind.funParams_;

% Automatic set up the starting points 
% For a Non-touching Cell!!!!
%[x, y] = findStartingPtOfWindowing(MD, iChan);
%windParams.StartPoint = [x(1), y(1)];

%windParams.StartPointPropag = false; % you might want this on for more normal cells I wanted it to remain constant
% Set the Windowing Parameters
windParams.SegProcessIndex = MD.getProcessIndex('MaskRefinementProcess');
windParams.ChannelIndex = iChan; %%%%%%%%%
windParams.ParaSize = windparasize;
windParams.PerpSize = windperpsize;
windParams.MinSize = 100; % min size of the object to window
windParams.StartContour = windstcontour;
%
windParams.MethodName = 'ConstantNumber';


% Note You can set the start window site - for instance I set it here to the
% neurite entrance point
%     SPIdx = veilStem(1).idxEnterNeurite;
%     [SPy,SPx] = ind2sub(MD.imSize_,SPIdx);
%     % windParams.StartPoint = [SPy SPx];
%     windParams.StartPoint = [SPx SPy]; %

parseProcessParams(MD.processes_{end}, windParams);

% Run the process
MD.processes_{end}.run();
MD.sanityCheck()


%% Run protrusion sampling process
protSamp = ProtrusionSamplingProcess(MD);
MD.addProcess(protSamp);
%protParams.OutputDirectory = [MD.outputDirectory_ filesep 'protrusion_samples_' outName];
% % Associate the threshold process to the package
%MD.packages_{2}.setProcess(3, protSamp);

protSampParams = protSamp.funParams_;
parseProcessParams(MD.processes_{end}, protSampParams);

% Run the process
MD.processes_{end}.run();
MD.sanityCheck()

%% WindowSamplingProcess
winSamp = WindowSamplingProcess(MD);
MD.addProcess(winSamp);

%protParams.OutputDirectory = [MD.outputDirectory_ filesep 'protrusion_samples_' outName];
% % Associate the threshold process to the package
%MD.packages_{2}.setProcess(4, winSamp);

winSampParams = winSamp.funParams_;
winSampParams.ChannelIndex = {windsamchindex}; %% [1;2;3;4];   %%%%%

parseProcessParams(MD.processes_{end}, winSampParams);

% Run the process
MD.processes_{end}.run();
MD.sanityCheck()

%% sampleMovieWindowsWithN
%tmp = sampleMovieWindowsWithN(MD);



%%
disp('== ML_DoExtSegRefProtWindWithN is done! ==')


end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5








    