% QZ change all 'numPerm', 100 to 'numPerm', 2 to make it running a little
% faster.
% QZ change all 'lagMax0', 50 to 'lagMax0', 5.
% QZ change all 'Arp3' & 'Apr3'(typo) to 'Actin' b/c wrong labelling.

%% example_pipeline_DX_FPAEME_2chan.m
%% The only input for this pipeline is the movieList object loaded to workspace.
% Jungsik Noh, 2018/05/08
% Qiongjing (Jenny) Zou,  JUN2018

%%  Diagnostic and Xcorrelation Analysis (DX)


%%  map Diagnostic plots and XCorrelation analysis for all windows
tic

maxLayer = 2

for i=1:numel(ML.movieDataFile_)
% movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName
    ML_DX_1chan(ML, i, maxLayer, 'Actin', 1, 'mapDescriptivesF', 'mapCrossCorrF')
    ML_DX_1chan(ML, i, maxLayer, 'mDia1', 2, 'mapDescriptivesF', 'mapCrossCorrF') % QZ diff from 1chan
    
    
    
%     ML_DX_xcorr2chan(ML, i, maxLayer, 'mDia1', 'Actin', 2, 1, 'mapCrossCorrF')  % QZ new for 2chan
    
%     function ML_DX_xcorr2chan(ML, MDindex, layerMax, chNametag1, chNametag2, ...
%         iChan1, iChan2, mapXCorrDirName)
    %  corr(iChan1_t+lag, iChan2_t)
    
    % QZ does it matter to put 'mDia1' or 'Actin' first???
    MDindex = i;
    layerMax = maxLayer;
    chNametag1 = 'mDia1';
    chNametag2 = 'Actin';
    iChan1 = 2;
    iChan2 = 1;
    mapXCorrDirName = 'mapCrossCorrF';
       
    % get movieData at the MDindex-th location
    load(fullfile(ML.movieDataFile_{MDindex}))      
    pause(1)                                             
    
    disp('===========================')
    disp(['======= MDindex = ', num2str(MDindex), ' ======='])
    disp('===========================')
    
    
    %%  outputDir
    figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName)
    
    
    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
    
    
    %% xcorr ch1, ch2
    
 
    chan1Name = chNametag1
    chan2Name = chNametag2

    
    mapXcorrCurvePermutation(MD, iChan1, iChan2, chan1Name, chan2Name,layerMax,figuresDir, ...
        'impute', 0, 'parpoolNum', 2, 'WithN', 0, 'numPerm', 2, ...
        'movingAvgSmoothing', 0, 'subFrames', subFr) %QZ I changed 'parpoolNum' from 12 to 2.
    
    
    for k=1:300; close(figure(k)); end
    
%     end

end

toc

%%  DX analysis summary at the ML level

maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,1,0,'Actin','Vel', maxLayer, ...
    'mapDescriptivesF', 'mapCrossCorrF', 'lagMax0', 5, ...
    'outDirName', 'Xcf_ch1ch0_F_lagMax50') % QZ 'lagMax0', 50 won't work for small # of frames, try # of frames/2
%%
maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,2,0,'mDia1','Vel', maxLayer, ...
    'mapDescriptivesF', 'mapCrossCorrF', 'lagMax0', 5, ...
    'outDirName', 'Xcf_ch2ch0_F_lagMax50')
%
maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,2,1,'mDia1','Actin', maxLayer, ...
    'mapDescriptivesF', 'mapCrossCorrF', 'lagMax0', 5, ...
    'outDirName', 'Xcf_ch2ch1_F_lagMax50')


%%  end: DX




%%  Fluctuation Profiling Around Edge Motion Events (FPAEME)


%%  ML_phaseDescriptives
tic
maxLayer = 2
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_phaseDescriptives(ML, i, maxLayer, 'Actin', 1, 'phaseDescriptives0p5mRL5')
end
toc
%  FPAME summary at the ML level

maxLayer = 2

MLsummary_FluctuationProfiling(ML, 'Actin', maxLayer, 'phaseDescriptives0p5mRL5', ...
    'phaseDesc_Arp3_0p5mRL5_lagMax20', 'lagMax0', 20)


%  ML_phaseDescriptives
tic
maxLayer = 2
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_phaseDescriptives(ML, i, maxLayer, 'mDia1', 2, 'phaseDescriptives0p5mRL5') % QZ diff, for chan No. 2.
end
toc
%  FPAME summary at the ML level

maxLayer = 2

MLsummary_FluctuationProfiling(ML, 'mDia1', maxLayer, 'phaseDescriptives0p5mRL5', ...
    'phaseDesc_mDia1_0p5mRL5_lagMax20', 'lagMax0', 20)



%%  end: FPAME




%%  ML_meanCV_alongDepth

tic
for i=1:numel(ML.movieDataFile_)
%     ML_meanCV_alongDepth(ML, i, 5, 'Actin', 1, 'mapDescriptives_alongDepth') % QZ Changed 'Apr3' (typo) to 'Actin'
    
%     function ML_meanCV_alongDepth(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
    % ML_meanCV_alongDepth_LB RUN mapDescriptives_meanCV_alongDepth() at ML
    % level.
    %
    % J Noh, 2018/01/30.
    
    MDindex = i;
    maxLayer = 5;
    chNametag = 'Actin';
    iChan = 1;
    outDirName = 'mapDescriptives_alongDepth';
    
    
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
    

    %%
    omitWin = []
    
    
    %%
    
    chanName = chNametag;
    
    mapDescriptives_meanCV_alongDepth(MD, iChan,maxLayer, chanName,figuresDir, ...
        'impute',0,'WithN',0, 'subFrames', subFr, 'omittedWindows', omitWin)
    
    for k=1:50; close(figure(k)); end
    
%     end
    
    
end
toc

%  MLsummary_alongDepth

MLsummary_alongDepth(ML,1,5,'Actin','mapDescriptives_alongDepth', ...
     'outDirName', 'topolayersMeanSD_Chan1')


%%
 
 

%%  for i=1:400; close(figure(i)); end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  movieData level parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
function ML_DX_1chan(ML, MDindex, maxLayer, chNametag, iChan, ...
    mapDescriptivesDirName, mapXCorrDirName)
% ML_DX_1chan RUN mapDescriptives_OneChan() and
% mapXcorrCurvePermutation_Vel() at the movieList level.
%
%
% Jungsik Noh, 2018/01/29.



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


%%  velocity Descriptives

iChan0 = 0;
chan0Name = 'Vel';
chan0Title = 'Velocity (nm/sec)';

mapDescriptives_OneChan(MD, iChan0, 1, chan0Name, chan0Title, figuresDir, 'adf', 1, ...
     'omittedWindows', [], 'subFrames', subFr, 'numPerm', 2, 'topograph', 'on', ...
     'parpoolNum', 2) %QZ I changed 'parpoolNum' from 6 to 2.



%%  Channel Descriptives

%iChan = 1;
chanName = chNametag;
chanTitle = chanName;
 
mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, 'adf', 1, ...
     'omittedWindows', [], 'subFrames', subFr, 'numPerm', 2, 'topograph', 'on', ...
     'parpoolNum', 2, 'WithN', 0) %QZ I changed 'parpoolNum' from 6 to 2.

for k=1:300; close(figure(k)); end           

            
%%  specify outputDir (cross correlation analysis)

%figuresDir = fullfile(MD.outputDirectory_, 'mapCrossCorrF');
figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);


%% xcorr chan vs velocity

layerMax = maxLayer;
%
mapXcorrCurvePermutation_Vel(MD, iChan, chanName, layerMax, figuresDir, 'impute', 0, ...
    'subFrames', subFr, 'parpoolNum', 2, 'WithN', 0, 'numPerm', 2) %QZ I changed 'parpoolNum' from 6 to 2.


for k=1:300; close(figure(k)); end

end


%%



%%
function ML_phaseDescriptives(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
% ML_phaseDescriptives_LB RUN phaseDescriptives_OneChan() at the ML level.
%
% Jungsik Noh, 2018/01/30.


%% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      
pause(1)    


%%  specify outputDir

%figuresDir = fullfile(MD.outputDirectory_, 'phaseDescriptives0p2mRL5_LB')
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
 
omitWin = []

%% phaseMasking to find a proper smoothing parameter

smParamTh = 0.5

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
end


%%


%%  EOF
