%% example_pipeline_DX_FPAEME_2chan_LB.m
%% The only input for this pipeline is the movieList object loaded to workspace.
% Jungsik Noh, 2018/05/08


%%  Diagnotics and Xcorrelation Analysis (DX)

%%  Identifying quiescent windows by using Ljung-Box test

for i=1:numel(ML.movieDataFile_)
    ML_Vel_classify_LB(ML, i, 'mapDescriptives_Vel_LB')
end



%%  map Diagnostic and XCorrelation analysis (DX) for only non-quiescent windows
tic

maxLayer = 2

for i=1:numel(ML.movieDataFile_)
% movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName
    ML_DX_1chan_LB(ML, i, maxLayer, 'Arp3', 1, 'mapDescriptives_LB', 'mapCrossCorr_LB')
    ML_DX_1chan_LB(ML, i, maxLayer, 'mDia1', 2, 'mapDescriptives_LB', 'mapCrossCorr_LB')
    ML_DX_xcorr2chan_LB(ML, i, maxLayer, 'mDia1', 'Arp3', 2, 1, 'mapCrossCorr_LB')
end

toc

%%  DX analysis summary at the ML level

maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,1,0,'Arp3','Vel', maxLayer, ...
    'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 40, ...
    'outDirName', 'Xcf_ch1ch0_LB_lagMax40')
%
maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,2,0,'mDia1','Vel', maxLayer, ...
    'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 40, ...
    'outDirName', 'Xcf_ch2ch0_LB_lagMax40')
%
maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,2,1,'mDia1','Arp3', maxLayer, ...
    'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 40, ...
    'outDirName', 'Xcf_ch2ch1_LB_lagMax40')


%%  end: DX




%%  Fluctuation Profiling Around Edge Motion Events (FPAEME)


%%  ML_phaseDescriptives
tic
maxLayer = 2
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_phaseDescriptives_LB(ML, i, maxLayer, 'Arp3', 1, 'phaseDescriptives0p5mRL5_LB')
end
toc
%  FPAME summary at the ML level

maxLayer = 2

MLsummary_FluctuationProfiling(ML, 'Arp3', maxLayer, 'phaseDescriptives0p5mRL5_LB', ...
    'phaseDesc_Arp3_LB_0p5mRL5_lagMax20', 'lagMax0', 20)


%  ML_phaseDescriptives
tic
maxLayer = 2
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_phaseDescriptives_LB(ML, i, maxLayer, 'mDia1', 2, 'phaseDescriptives0p5mRL5_LB')
end
toc
%  FPAME summary at the ML level

maxLayer = 2

MLsummary_FluctuationProfiling(ML, 'mDia1', maxLayer, 'phaseDescriptives0p5mRL5_LB', ...
    'phaseDesc_mDia1_LB_0p5mRL5_lagMax20', 'lagMax0', 20)



%%  end: FPAME




%%  ML_meanCV_alongDepth_LB

tic
for i=1:numel(ML.movieDataFile_)
    ML_meanCV_alongDepth_LB(ML, i, 5, 'Apr3', 1, 'mapDescriptives_alongDepth_LB')
end
toc

%  MLsummary_alongDepth

MLsummary_alongDepth(ML,1,5,'Arp3','mapDescriptives_alongDepth_LB', ...
     'outDirName', 'topolayersMeanSD_LB_Chan1')


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

%figuresDir = fullfile(MD.outputDirectory_, 'mapDescriptivesF');
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
     'omittedWindows', [], 'subFrames', subFr, 'numPerm', 100, 'topograph', 'on', ...
     'parpoolNum', 6)



%%  Channel Descriptives

%iChan = 1;
chanName = chNametag;
chanTitle = chanName;
 
mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, 'adf', 1, ...
     'omittedWindows', [], 'subFrames', subFr, 'numPerm', 100, 'topograph', 'on', ...
     'parpoolNum', 6, 'WithN', 0)

for k=1:300; close(figure(k)); end           

            
%%  specify outputDir (cross correlation analysis)

%figuresDir = fullfile(MD.outputDirectory_, 'mapCrossCorrF');
figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);


%% xcorr chan vs velocity

layerMax = maxLayer;
%
mapXcorrCurvePermutation_Vel(MD, iChan, chanName, layerMax, figuresDir, 'impute', 0, ...
    'subFrames', subFr, 'parpoolNum', 6, 'WithN', 0, 'numPerm', 100)


for k=1:300; close(figure(k)); end

end


%%
function ML_DX_1chan_LB(ML, MDindex, maxLayer, chNametag, iChan, ...
    mapDescriptivesDirName, mapXCorrDirName)
% ML_DX_1chan_LB RUN mapDescriptives_OneChan() and
% mapXcorrCurvePermutation_Vel() at the movieList level after excluding
% quiescent windows.
%
%
% Jungsik Noh, 2018/01/29.



%% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      
pause(1)                                        
                                                

%%  specify outputDir

%figuresDir = fullfile(MD.outputDirectory_, 'mapDescriptives_LB');
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

%figuresDir = fullfile(MD.outputDirectory_, 'mapCrossCorr_LB');
figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);


%% xcorr chan vs velocity

layerMax = maxLayer;
%
mapXcorrCurvePermutation_Vel(MD, iChan, chanName, layerMax, figuresDir, 'impute', 0, ...
    'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr, ...
    'WithN', 0, 'numPerm', 100, 'lagMax', 40)


for k=1:300; close(figure(k)); end

end



%%
function ML_DX_xcorr2chan(ML, MDindex, layerMax, chNametag1, chNametag2, ...
    iChan1, iChan2, mapXCorrDirName)
%  corr(iChan1_t+lag, iChan2_t)

% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      %ML.getMovies();
pause(1)                                                %MD = ML.getMovie(MDindex);

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

%iChan1 = 11
%iChan2 = 2
chan1Name = chNametag1
chan2Name = chNametag2
%layerMax = 2;

mapXcorrCurvePermutation(MD, iChan1, iChan2, chan1Name, chan2Name,layerMax,figuresDir, ...
   'impute', 0, 'parpoolNum', 12, 'WithN', 0, 'numPerm', 100, ...
   'movingAvgSmoothing', 0, 'subFrames', subFr)
    

for k=1:300; close(figure(k)); end

end


%%
function ML_DX_xcorr2chan_LB(ML, MDindex, layerMax, chNametag1, chNametag2, ...
    iChan1, iChan2, mapXCorrDirName)
%  corr(iChan1_t+lag, iChan2_t)

% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      %ML.getMovies();
pause(1)                                                %MD = ML.getMovie(MDindex);

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

%%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
%%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.

velAnalName = 'mapDescriptives_Vel_LB';
indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
raw = load(indPath);

omitWin = find(raw.indActive == 0);
disp(omitWin)



%% xcorr ch1, ch2

%iChan1 = 11
%iChan2 = 2
chan1Name = chNametag1
chan2Name = chNametag2
%layerMax = 2;

mapXcorrCurvePermutation(MD, iChan1, iChan2, chan1Name, chan2Name,layerMax,figuresDir, ...
   'impute', 0, 'parpoolNum', 12, 'WithN', 0, 'numPerm', 100, 'omittedWindows', omitWin, ...
   'movingAvgSmoothing', 0, 'subFrames', subFr)
    

for k=1:300; close(figure(k)); end

end



%%
function ML_Vel_classify_LB(ML, MDindex, outDirName)
% ML_Vel_classify_LB RUN mapDescriptives_Vel_LB.m at the movieList level.
%
%
% Jungsik Noh, 2018/01/29.


%%  get movieData
load(fullfile(ML.movieDataFile_{MDindex}));
pause(1)                                                
                                                

%%  specify outputDir

%figuresDir = fullfile(MD.outputDirectory_, 'mapDescriptives_Vel_LB');
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

mapDescriptives_Vel_LB(MD, figuresDir, 'impute',1,'movingAvgSmoothing',1, ...
    'figFlag','off', 'omittedWindows', [], 'subFrames', subFr, 'topograph', 'on')


for k=1:300; close(figure(k)); end

end


%%
function ML_phaseDescriptives_LB(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
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
end



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
function ML_meanCV_alongDepth_LB(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
% ML_meanCV_alongDepth_LB RUN mapDescriptives_meanCV_alongDepth() at ML
% level.
% 
% J Noh, 2018/01/30.

%% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      %ML.getMovies();
pause(1)                                                %MD = ML.getMovie(MDindex);

%%  outputDir
%figuresDir = fullfile(MD.outputDirectory_, 'mapDescriptives_alongDepth')
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


velAnalName = 'mapDescriptives_Vel_LB'
indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
raw = load(indPath);

omitWin = find(raw.indActive == 0)

%%  

chanName = chNametag;

mapDescriptives_meanCV_alongDepth(MD, iChan,maxLayer, chanName,figuresDir, ...
    'impute',0,'WithN', 0, 'subFrames', subFr, 'omittedWindows', omitWin)

for k=1:50; close(figure(k)); end

end



%%
function ML_meanCV_alongDepth(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
% ML_meanCV_alongDepth_LB RUN mapDescriptives_meanCV_alongDepth() at ML
% level.
% 
% J Noh, 2018/01/30.

%% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      %ML.getMovies();
pause(1)                                                %MD = ML.getMovie(MDindex);

%%  outputDir
%figuresDir = fullfile(MD.outputDirectory_, 'mapDescriptives_alongDepth')
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


%%  

chanName = chNametag;

mapDescriptives_meanCV_alongDepth(MD, iChan,maxLayer, chanName,figuresDir, ...
    'impute',0,'WithN',0, 'subFrames', subFr, 'omittedWindows', omitWin)

for k=1:50; close(figure(k)); end

end


 



%%  EOF
