% example_pipeline_DX_FPAEME_1chan_LB_actin.m
% Use MovieData cell array as input for all functions.
% Jungsik Noh, 2018/05/07
% Qiongjing (Jenny) Zou,  JUN2018

%% creat MDs from raw images:
% 25 frames for each MD, 2018-07-02
path1 = 'C:\Users\Qiongjing Zou\Documents\MATLAB\example2movies\QZdata1';
fullpath1 = fullfile(path1,'images_actin');
c1 = Channel(fullpath1);
MD1 = MovieData(c1,path1); % path1 here is MD's outputDirectory_
MD1.setPath(path1); % set movieDataPath_, where to save .mat file

path2 = 'C:\Users\Qiongjing Zou\Documents\MATLAB\example2movies\QZdata2';
fullpath2 = fullfile(path2,'images_actin');
c2 = Channel(fullpath2);
MD2 = MovieData(c2,path2);
MD2.setPath(path2);

MDs=cell(1,2);
MDs{1} = MD1;
MDs{2} = MD2;
%% setup MDs, add 6 processes
for i = 1:2
    MDs{i}.setFilename(sprintf('moviedata%dactin.mat',i));
    MDs{i}.pixelSize_ = 108;
    MDs{i}.timeInterval_ = 5;
    MDs{i}.sanityCheck; % add nFrames_, imSize_, and reader to the object. also save movie by obj.save()
    
    % add 6 processes.
    MDs{i}.reset()
    process = MultiScaleAutoSegmentationProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run() % did obj.getOwner().save()
    
    process = MaskRefinementProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    
    process = ProtrusionProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    % window will appear to ask to choose 'Mask Refinement' as mask process.
    
    process = WindowingProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    % window will appear to ask to choose 'Mask Refinement' as mask process.
    
    process = ProtrusionSamplingProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    
    process = WindowSamplingProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
end

%% no need to run above script every time, can load from saved mat files:
MDs=cell(1,2);
folder='C:\Users\Qiongjing Zou\Documents\MATLAB\example2movies';
for i = 1:2
load([folder filesep sprintf('QZdata%d', i) filesep sprintf('movieData%dactin.mat',i)]);
MDs{i}=MD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Diagnotics and Xcorrelation Analysis (DX):

%%  Identifying quiescent windows by using Ljung-Box test!
% QZ this section is new for LB pipeline.
% QZ the only useful thing from this LB test is the indActive_windowIndex.mat in its outputdirectory

for i=1:numel(MDs)  
    %%  REPEATED Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    %%  run mapDescriptives_Vel_LB with options
    %  specify outputDir
    outDirName = 'mapDescriptives_Vel_LB';
    figuresDir = fullfile(MDs{i}.outputDirectory_, outDirName);
    
    %  mapDescriptives_Vel_LB(MD, figuresDir, varargin)
    mapDescriptives_Vel_LB(MDs{i}, figuresDir, 'impute',1,'movingAvgSmoothing',1, ...
        'figFlag','off', 'omittedWindows', [], 'subFrames', subFr, 'topograph', 'on')

    close all
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MapDDX:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  map Diagnostic and XCorrelation analysis (DX) for only non-quiescent windows
tic

maxLayer = 2

for i=1:numel(MDs)
% movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName    
%     ML_DX_1chan_LB(ML, i, maxLayer, 'Actin', 1, 'mapDescriptives_LB', 'mapCrossCorr_LB')
    
%     function ML_DX_1chan_LB(ML, MDindex, maxLayer, chNametag, iChan, ...
%         mapDescriptivesDirName, mapXCorrDirName)
    % ML_DX_1chan_LB RUN mapDescriptives_OneChan() and
    % mapXcorrCurvePermutation_Vel() at the movieList level after excluding
    % quiescent windows.

    chNametag = 'Actin';
    iChan = 1;
       
    %%  REPEATED Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    %%  REPEATED exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    % QZ this section is new in LB pipeline.
    
    velAnalName = 'mapDescriptives_Vel_LB';
    indPath = fullfile(MDs{i}.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath); % QZ load a .mat file calculated from LB test section
    
    omitWin = find(raw.indActive == 0);
    disp(omitWin)
    
    
    %%  specify outputDir
    mapDescriptivesDirName = 'mapDescriptives_LB';
    figuresDir = fullfile(MDs{i}.outputDirectory_, mapDescriptivesDirName);
    
    %%  velocity Descriptives!
    
    iChan0 = 0;
    chan0Name = 'Vel';
    chan0Title = 'Velocity (nm/sec)';
    
    mapDescriptives_OneChan(MDs{i}, iChan0, 1, chan0Name, chan0Title, figuresDir, 'adf',1,'impute',0, ...
        'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr,'numPerm', 2, 'topograph', 'on') % QZ change all 'numPerm', 100 to 'numPerm', 2 to make it running a little faster.
    
    
    
    %%  Chan Descriptives!

    chanName = chNametag;
    chanTitle = chanName;
    
    mapDescriptives_OneChan(MDs{i}, iChan, maxLayer, chanName, chanTitle, figuresDir, 'adf', 1, ...
        'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr,'numPerm', 2, 'topograph', 'on', 'WithN', 0)
    
    close all    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  specify outputDir of correlation analysis
    mapXCorrDirName = 'mapCrossCorr_LB';      
    figuresDir = fullfile(MDs{i}.outputDirectory_, mapXCorrDirName);
        
    %% xcorr chan vs velocity!
    
    layerMax = maxLayer;
    
    mapXcorrCurvePermutation_Vel(MDs{i}, iChan, chanName, layerMax, figuresDir, 'impute', 0, ...
        'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr, ...
        'WithN', 0, 'numPerm', 2, 'lagMax', 4) % QZ change 'lagMax' from 40
    
    close all
    
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DX analysis summary at the ML level

maxLayer = 2

MDsSummary_XcorrCurvesVelAcf(MDs,1,0,'Actin','Vel', maxLayer, ...
    'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 4, ...
    'outDirName', 'SummaryXcf_ch1ch0_LB_lagMax4') % QZ change 'lagMax0' from 40
% old: 'outDirName', 'Xcf_ch1ch0_LB_lagMax40'


%% end: DX 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FPEAME:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Fluctuation Profiling Around Edge Motion Events (FPEAME)


%%  ML_phaseDescriptives_LB
tic
maxLayer = 2
for i=1:numel(MDs)
    % check smParam and minimumRunLength within ML_ function.
%     ML_phaseDescriptives_LB(ML, i, maxLayer, 'Actin', 1, 'phaseDescriptives0p4mRL5_LB')
    
%     function ML_phaseDescriptives_LB(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
    % ML_phaseDescriptives_LB RUN phaseDescriptives_OneChan() at the ML level.
    %

    MDindex = i;
    chNametag = 'Actin';
    iChan = 1;    
    
    %%  REPEATED Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
   
    %%  RPEATED exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    % QZ this section is new in LB pipeline.
    
    velAnalName = 'mapDescriptives_Vel_LB';
    indPath = fullfile(MDs{i}.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath); % QZ load a .mat file calculated from LB test section
    
    omitWin = find(raw.indActive == 0)    
    
    %%  specify outputDir
    outDirName = 'phaseDescriptives_LB';
    figuresDir = fullfile(MDs{i}.outputDirectory_, outDirName);
    
    %% phaseMasking to find a proper smoothing parameter
    
    smParamTh = 0.4
    
    [protMask1, retMask1] = phaseMasking(MDs{i}, smParamTh, figuresDir, 'minimumRunLength', 5, ...
        'subFrames', subFr);
    
    %%  Retraction / specify sampling bandwidth
    Mask0 = retMask1;
    samplingBw = 20
    
    %%  Chan Descriptives
    
    chanName = [chNametag, '-ret'];
    chanTitle = chanName;

    phaseDescriptives_OneChan(MDs{i}, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0) 
    
    %%  Protrusion
    Mask0 = protMask1;
    
    %%  Chan Descriptives
    
    chanName = [chNametag, '-prot'];
    chanTitle = chanName;

    phaseDescriptives_OneChan(MDs{i}, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    %%  phase MaxMinVel
    
    chanName = chNametag;
    chanTitle = chanName;
 
    phaseDescriptives_MaxMinVel_OneChan(MDs{i}, iChan, maxLayer, chanName, chanTitle, ...
        smParamTh, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    close all

end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  FPAME summary at the ML level

maxLayer = 2

MDsSummary_FluctuationProfiling(MDs, 'Actin', maxLayer, 'phaseDescriptives_LB', ...
    'SummaryphaseDesc_0p4mRL5_LB_lagMax4', 'lagMax0', 4) % QZ change 'lagMax0' from 10
% old 'phaseDesc_0p4mRL5_LB_lagMax10'

close all

%%  end: FPAME

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%alongDepth:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ML_meanCV_alongDepth

tic
for i=1:numel(MDs)
%     ML_meanCV_alongDepth_LB(ML, i, 7, 'Actin', 1, 'mapDescriptives_alongDepth_LB')
    
    
%     function ML_meanCV_alongDepth_LB(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
    % ML_meanCV_alongDepth_LB RUN mapDescriptives_meanCV_alongDepth() at ML
    % level.
    %
    
    MDindex = i;
    maxLayer = 7;
    chNametag = 'Actin';
    iChan = 1;
    outDirName = 'mapDescriptives_alongDepth_LB';                                              
    
    %%  outputDir
    
    figuresDir = fullfile(MDs{i}.outputDirectory_, outDirName);
    
    %%  REPEATED Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MDs{i}.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end

    %%  RPEATED exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    % QZ this section is new in LB pipeline.
    
    velAnalName = 'mapDescriptives_Vel_LB';
    indPath = fullfile(MDs{i}.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath); % QZ load a .mat file calculated from LB test section
    
    omitWin = find(raw.indActive == 0)    
    
    %%
    
    chanName = chNametag;
    
    mapDescriptives_meanCV_alongDepth(MDs{i}, iChan,maxLayer, chanName,figuresDir, ...
        'impute',0,'WithN', 0, 'subFrames', subFr, 'omittedWindows', omitWin)
    
    close all
    
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MLsummary_alongDepth

MDsSummary_alongDepth(MDs,1,7,'Actin','mapDescriptives_alongDepth_LB', ...
     'outDirName', 'SummarytopolayersMeanSD_LB_Chan1')
 % old 'topolayersMeanSD_LB_Chan1'

close all
%%  EOF
