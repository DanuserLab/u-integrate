%% example_pipeline_DX_FPAEME_1chan_LB_actin.m
%% The only input for this pipeline is the movieList object loaded to workspace.
% Jungsik Noh, 2018/05/07
% Qiongjing (Jenny) Zou,  JUN2018


%%  Diagnotics and Xcorrelation Analysis (DX)

%%  Identifying quiescent windows by using Ljung-Box test
% QZ this section is new for LB pipeline.

for i=1:numel(ML.movieDataFile_)
%     ML_Vel_classify_LB(ML, i, 'mapDescriptives_Vel_LB')
    
    % function ML_Vel_classify_LB(ML, MDindex, outDirName)
    % ML_Vel_classify_LB RUN mapDescriptives_Vel_LB.m at the movieList level.
    %
    %
    % Jungsik Noh, 2018/01/29.
    
    MDindex = i;
    outDirName = 'mapDescriptives_Vel_LB';
    
    %%  get movieData
    load(fullfile(ML.movieDataFile_{MDindex}));
    pause(1)
    
    
    %%  specify outputDir

    figuresDir = fullfile(MD.outputDirectory_, outDirName);
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    %%  run mapDescriptives_Vel_LB with options
    
    %  mapDescriptives_Vel_LB(MD, figuresDir, varargin)
    mapDescriptives_Vel_LB(MD, figuresDir, 'impute',1,'movingAvgSmoothing',1, ...
        'figFlag','off', 'omittedWindows', [], 'subFrames', subFr, 'topograph', 'on')
    
    
    for k=1:300; close(figure(k)); end
    
    % end
    
end


%%  map Diagnostic and XCorrelation analysis (DX) for only non-quiescent windows
tic

maxLayer = 2

for i=1:numel(ML.movieDataFile_)
% movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName    
%     ML_DX_1chan_LB(ML, i, maxLayer, 'Actin', 1, 'mapDescriptives_LB', 'mapCrossCorr_LB')
    
%     function ML_DX_1chan_LB(ML, MDindex, maxLayer, chNametag, iChan, ...
%         mapDescriptivesDirName, mapXCorrDirName)
    % ML_DX_1chan_LB RUN mapDescriptives_OneChan() and
    % mapXcorrCurvePermutation_Vel() at the movieList level after excluding
    % quiescent windows.
    %
    %
    % Jungsik Noh, 2018/01/29.
    
    MDindex = i;
    chNametag = 'Actin';
    iChan = 1;
    mapDescriptivesDirName = 'mapDescriptives_LB';
    mapXCorrDirName = 'mapCrossCorr_LB';
    
    
    %% get movieData at the MDindex-th location
    load(fullfile(ML.movieDataFile_{MDindex}))
    pause(1)
    
    
    %%  specify outputDir
    
    figuresDir = fullfile(MD.outputDirectory_, mapDescriptivesDirName);
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    % QZ this section is new in LB pipeline.
    
    velAnalName = 'mapDescriptives_Vel_LB';
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0);
    disp(omitWin)
    
    
    %%  velocity Descriptives
    
    iChan0 = 0;
    chan0Name = 'Vel';
    chan0Title = 'Velocity (nm/sec)';
    
    mapDescriptives_OneChan(MD, iChan0, 1, chan0Name, chan0Title, figuresDir, 'adf',1,'impute',0, ...
        'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr,'numPerm', 100, 'topograph', 'on')
    
    
    
    %%  Chan Descriptives
    
    %iChan = 1;
    chanName = chNametag;
    chanTitle = chanName;
    
    mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, 'adf', 1, ...
        'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr,'numPerm', 100, 'topograph', 'on', 'WithN', 0)
    
    for k=1:300; close(figure(k)); end
    
    
    %%  specify outputDir of correlation analysis
    
    figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);
    
    
    %% xcorr chan vs velocity
    
    layerMax = maxLayer;
    %
    mapXcorrCurvePermutation_Vel(MD, iChan, chanName, layerMax, figuresDir, 'impute', 0, ...
        'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr, ...
        'WithN', 0, 'numPerm', 100, 'lagMax', 40)
    
    
    for k=1:300; close(figure(k)); end
    
%     end
    
end

toc

%%  DX analysis summary at the ML level

maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,1,0,'Actin','Vel', maxLayer, ...
    'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 40, ...
    'outDirName', 'Xcf_ch1ch0_LB_lagMax40')



%% end: DX 



%%  Fluctuation Profiling Around Edge Motion Events (FPEAME)


%%  ML_phaseDescriptives_LB
tic
maxLayer = 2
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
%     ML_phaseDescriptives_LB(ML, i, maxLayer, 'Actin', 1, 'phaseDescriptives0p4mRL5_LB')
    
%     function ML_phaseDescriptives_LB(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
    % ML_phaseDescriptives_LB RUN phaseDescriptives_OneChan() at the ML level.
    %
    % Jungsik Noh, 2018/01/30.
    
    MDindex = i;
    chNametag = 'Actin';
    iChan = 1;
    outDirName = 'phaseDescriptives0p4mRL5_LB';
    
    
    %% get movieData at the MDindex-th location
    load(fullfile(ML.movieDataFile_{MDindex}))
    pause(1)
    
    
    %%  specify outputDir
    
    figuresDir = fullfile(MD.outputDirectory_, outDirName);
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    
    %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    % QZ this section is new in LB pipeline.
    
    velAnalName = 'mapDescriptives_Vel_LB'
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0)
    
    
    %% phaseMasking to find a proper smoothing parameter
    
    smParamTh = 0.4
    
    [protMask1, retMask1] = phaseMasking(MD, smParamTh, figuresDir, 'minimumRunLength', 5, ...
        'subFrames', subFr);
    
    
    
    %%  Retraction / specify sampling bandwidth
    Mask0 = retMask1;
    samplingBw = 20
    
    
    %%  Chan Descriptives
    
    chanName = [chNametag, '-ret'];
    chanTitle = chanName;
    
    %
    phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    
    
    %%  Protrusion
    Mask0 = protMask1;
    
    
    %%  Chan Descriptives
    
    chanName = [chNametag, '-prot'];
    chanTitle = chanName;
    
    %
    phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    
    
    %%  phase MaxMinVel
    
    chanName = chNametag;
    chanTitle = chanName;
    
    
    phaseDescriptives_MaxMinVel_OneChan(MD, iChan, maxLayer, chanName, chanTitle, ...
        smParamTh, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    
    
    for k=1:300; close(figure(k)); end
    
    
    %%
    %     end

end
toc


%%  FPAME summary at the ML level

maxLayer = 2

MLsummary_FluctuationProfiling(ML, 'Actin', maxLayer, 'phaseDescriptives0p4mRL5_LB', ...
    'phaseDesc_0p4mRL5_LB_lagMax10', 'lagMax0', 10)


%%  end: FPAME



%% ML_meanCV_alongDepth

tic
for i=1:numel(ML.movieDataFile_)
%     ML_meanCV_alongDepth_LB(ML, i, 7, 'Actin', 1, 'mapDescriptives_alongDepth_LB')
    
    
%     function ML_meanCV_alongDepth_LB(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
    % ML_meanCV_alongDepth_LB RUN mapDescriptives_meanCV_alongDepth() at ML
    % level.
    %
    % J Noh, 2018/01/30.
    
    MDindex = i;
    maxLayer = 7;
    chNametag = 'Actin';
    iChan = 1;
    outDirName = 'mapDescriptives_alongDepth_LB';
    
    %% get movieData at the MDindex-th location
    load(fullfile(ML.movieDataFile_{MDindex}))    
    pause(1)                                               
    
    %%  outputDir
    
    figuresDir = fullfile(MD.outputDirectory_, outDirName);
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    % QZ this section is new in LB pipeline.
    
    velAnalName = 'mapDescriptives_Vel_LB'
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0)
    
    %%
    
    chanName = chNametag;
    
    mapDescriptives_meanCV_alongDepth(MD, iChan,maxLayer, chanName,figuresDir, ...
        'impute',0,'WithN', 0, 'subFrames', subFr, 'omittedWindows', omitWin)
    
    for k=1:50; close(figure(k)); end
    
%     end
    
end
toc


%%  MLsummary_alongDepth

MLsummary_alongDepth(ML,1,7,'Actin','mapDescriptives_alongDepth_LB', ...
     'outDirName', 'topolayersMeanSD_LB_Chan1')
 
 

%% 
%% for i=1:100; close(figure(i)); end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  movieData level parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%




%%



%%



%%

 

%%  EOF
