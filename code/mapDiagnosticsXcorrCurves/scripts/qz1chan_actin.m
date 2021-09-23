%% example_pipeline_DX_FPAEME_1chan_actin.m
%% The only input for this pipeline is the movieList object loaded to workspace.
% Jungsik Noh, 2018/05/07
% Qiongjing (Jenny) Zou,  JUN2018


%%  Diagnotics and Xcorrelation Analysis (DX)

%% map Diagnostic plots and XCorrelation analysis for all windows
tic

maxLayer = 2

for i=1:numel(ML.movieDataFile_) % QZ there are 2 movieData.mat files, each has 2 channels.
% movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName
%     ML_DX_1chan(ML, i, maxLayer, 'Actin', 1, 'mapDescriptivesF', 'mapCrossCorrF')
    
%     function ML_DX_1chan(ML, MDindex, maxLayer, chNametag, iChan, ...
%         mapDescriptivesDirName, mapXCorrDirName)

    % ML_DX_1chan RUN mapDescriptives_OneChan() and
    % mapXcorrCurvePermutation_Vel() at the movieList level.
    %
    %
    % Jungsik Noh, 2018/01/29.
    
    MDindex = i;
    chNametag = 'Actin';
  
%     mapDescriptivesDirName = 'mapDescriptivesF';
%     mapXCorrDirName = 'mapCrossCorrF';
    mapDescriptivesDirName = 'S1F1mapDesc'; % QZ S1 is script 1: 1chan_actin, F1 is folder1
    mapXCorrDirName = 'S1F2mapXcorr';    
    
    
    %% get movieData at the MDindex-th location
    load(fullfile(ML.movieDataFile_{MDindex})) % QZ load MD
    pause(1)
    
    
    %%  specify outputDir
    
    figuresDir = fullfile(MD.outputDirectory_, mapDescriptivesDirName);
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    % QZ what is purpose of this section, I did not see a subFr.mat???
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    %%  velocity Descriptives
    
    iChan0 = 0;
    chan0Name = 'Vel';
    chan0Title = 'Velocity (nm/sec)';
    % QZ for Velocity, maxLayer is 1.
     
    % mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, varargin)
    mapDescriptives_OneChan(MD, iChan0, 1, chan0Name, chan0Title, figuresDir, 'adf', 1, ...
        'omittedWindows', [], 'subFrames', subFr, 'numPerm', 100, 'topograph', 'on', 'parpoolNum', 2) %QZ I changed 'parpoolNum' from 6 to 2.
   
    
    
    
    %%  Channel Descriptives
    
    iChan = 1;
    chanName = chNametag;
    chanTitle = chanName;
    
    mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, 'adf', 1, ...
        'omittedWindows', [], 'subFrames', subFr, 'numPerm', 100, 'topograph', 'on', 'parpoolNum', 2, 'WithN', 0) %QZ I changed 'parpoolNum' from 6 to 2.
    
    
    % for k=1:300; close(figure(k)); end
    close all % QZ
    
    %%  specify outputDir (cross correlation analysis)
    
    figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);
    
    
    %% xcorr chan vs velocity
    
    % mapXcorrCurvePermutation_Vel(MD, iChan1, chan1Name, layerMax, figuresDir, varargin) 
    mapXcorrCurvePermutation_Vel(MD, iChan, chanName, maxLayer, figuresDir, 'impute', 0, ...
        'subFrames', subFr, 'parpoolNum', 2, 'WithN', 0, 'numPerm', 100, 'lagMax', 40) %QZ I changed 'parpoolNum' from 6 to 2.
    
    
     % for k=1:300; close(figure(k)); end
    close all % QZ
    
    % end
    
end

toc


%% DX analysis summary at the ML level

maxLayer = 2

% MLsummary_XcorrCurvesVelAcf(ML, iChan1, iChan2, chan1Name, chan2Name, maxLayer, ...
%    analNameAcf, analNameXcf, varargin)
% MLsummary_XcorrCurvesVelAcf uses the output figures from
% mapDescriptives_OneChan and mapXcorrCurvePermutation_Vel as input ('mapDescriptivesF', 'mapCrossCorrF')

% MLsummary_XcorrCurvesVelAcf(ML,1,0,'Actin','Vel', maxLayer, ...
%     'mapDescriptivesF', 'mapCrossCorrF', 'lagMax0', 4, ...
%     'outDirName', 'Xcf_ch1ch0_F_lagMax40') % QZ change 'lagMax0', 40 to 'lagMax0', 4.; Velocity (iChan is 0) as chan2.
MLsummary_XcorrCurvesVelAcf(ML,1,0,'Actin','Vel', maxLayer, ...
    mapDescriptivesDirName, mapXCorrDirName, 'lagMax0', 4, ...
    'outDirName', 'S1F3mapDDX_summary') % QZ change 'lagMax0', 40 to 'lagMax0', 4.; Velocity (iChan is 0) as chan2.
  
%%
close all % QZ added

%% end: DX 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Fluctuation Profiling Around Edge Motion Events (FPEAME)


%%  ML_phaseDescriptives
% QZ was ML_phaseDescriptives_LB, I feel it is wrong
tic
maxLayer = 2
for i=1:numel(ML.movieDataFile_)
    % Check the tuning parameters of smParam, minimumRunLength, samplingBw within ML_phaseDescriptives function.
%     ML_phaseDescriptives(ML, i, maxLayer, 'Actin', 1, 'phaseDescriptives0p4mRL5')
    
%     function ML_phaseDescriptives(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
    % ML_phaseDescriptives RUN phaseDescriptives_OneChan() at the ML level.
    %
    % Jungsik Noh, 2018/01/30.
    
    MDindex = i;
    chNametag = 'Actin';
    iChan = 1;
%     outDirName = 'phaseDescriptives0p4mRL5';
    outDirName = 'S1F4phaseDesc';
    
    %% get movieData at the MDindex-th location
    load(fullfile(ML.movieDataFile_{MDindex}))
    pause(1)
    
    
    %%  specify outputDir
    
    figuresDir = fullfile(MD.outputDirectory_, outDirName);
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    % QZ this part was repeated.
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    %%  omitWin
    
    omitWin = [];
    
    %% phaseMasking to find a proper smoothing parameter
    
    smParamTh = 0.4
    
    % [protMask1, retMask1] = phaseMasking(MD, smParamTh, figuresDir, varargin)
    [protMask1, retMask1] = phaseMasking(MD, smParamTh, figuresDir, 'minimumRunLength', 5, ...
        'subFrames', subFr); % QZ MD was loaded to workspace.
    
    
    
    %%  Retraction / specify sampling bandwidth
    Mask0 = retMask1;
    samplingBw = 20
    
    
    %%  Chan Descriptives
    
    chanName = [chNametag, '-ret']; % QZ 'ret' is retraction
    chanTitle = chanName;
    
    %
    % phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask, samplingBw, figuresDir, varargin)
    phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    
    
    %%  Protrusion
    Mask0 = protMask1;
    
    
    %%  Chan Descriptives
    
    chanName = [chNametag, '-prot']; % QZ 'prot' is protrusion
    chanTitle = chanName;
    
    %
    % phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask, samplingBw, figuresDir, varargin)
    phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    
    
    %%  phase MaxMinVel
    
    chanName = chNametag;
    chanTitle = chanName;
    
    % phaseDescriptives_MaxMinVel_OneChan(MD, iChan, maxLayer, chanName, chanTitle, ...
    %        smParamTh, samplingBw, figuresDir, varargin)
    phaseDescriptives_MaxMinVel_OneChan(MD, iChan, maxLayer, chanName, chanTitle, ...
        smParamTh, samplingBw, figuresDir, ...
        'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)
    
    
    
    % for k=1:300; close(figure(k)); end
    close all % QZ
      
    % end

end
toc


%% FPAME summary at the ML level
maxLayer = 2

% MLsummary_FluctuationProfiling(ML, chanName, maxLayer, analNameDesc, outDirName, varargin)
% MLsummary_FluctuationProfiling(ML, 'Actin', maxLayer, 'phaseDescriptives0p4mRL5', ...
%     'phaseDesc3_0p4mRL5_lagMax10', 'lagMax0', 10)
MLsummary_FluctuationProfiling(ML, 'Actin', maxLayer, 'S1F4phaseDesc', ...
    'S1F5FPAEME_summary', 'lagMax0', 10)

%%
close all % QZ added

%%  end: FPAME

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ML_meanCV_alongDepth

tic
for i=1:numel(ML.movieDataFile_)
%     ML_meanCV_alongDepth(ML, i, 7, 'Actin', 1, 'mapDescriptives_alongDepth')
    
    % function ML_meanCV_alongDepth(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
    % ML_meanCV_alongDepth RUN mapDescriptives_meanCV_alongDepth() at ML
    % level.
    %
    % J Noh, 2018/01/30.
    
    MDindex = i;
    maxLayer = 7;
    chNametag = 'Actin';
    iChan = 1;
%     outDirName = 'mapDescriptives_alongDepth';
    outDirName = 'S1F6mapDesc_alongDepth';
       
    %% get movieData at the MDindex-th location
    load(fullfile(ML.movieDataFile_{MDindex}))      
    pause(1)                                              
    
    
    %%  outputDir
    
    figuresDir = fullfile(MD.outputDirectory_, outDirName);
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    % QZ this part was repeated 3rd time.
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    
    %%
    
    chanName = chNametag;
    
    % mapDescriptives_meanCV_alongDepth(MD, iChan, maxLayer, chanName, figuresDir, varargin)
    mapDescriptives_meanCV_alongDepth(MD, iChan, maxLayer, chanName, figuresDir, ...
        'impute',0,'WithN',0, 'subFrames', subFr, 'omittedWindows', [])
    
    % for k=1:50; close(figure(k)); end
    close all % QZ
    
    % end
    
    
end
toc


%%  MLsummary_alongDepth

% MLsummary_alongDepth(ML, iChan, maxLayer, chanName, analFolderName, ...
%        varargin)
% MLsummary_alongDepth(ML,1,7,'Actin','mapDescriptives_alongDepth', ...
%      'outDirName', 'topolayersMeanSD3_Chan1')
MLsummary_alongDepth(ML,1,7,'Actin','S1F6mapDesc_alongDepth', ...
     'outDirName', 'S1F7alongDepth_summary')
 
%%
close all % QZ added

%% 
%% for i=1:100; close(figure(i)); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  movieData level parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%



%%



%%



 

%%  EOF
