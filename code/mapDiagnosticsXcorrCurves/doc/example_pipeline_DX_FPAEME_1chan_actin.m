%% example_pipeline_DX_FPAEME_1chan_actin.m
%% The only input for this pipeline is the movieList object loaded to workspace.
% Jungsik Noh, 2018/05/07


%%  Diagnotics and Xcorrelation Analysis (DX)

%% map Diagnostic plots and XCorrelation analysis for all windows
tic

maxLayer = 2

for i=1:numel(ML.movieDataFile_)
% movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName
    ML_DX_1chan(ML, i, maxLayer, 'Actin', 1, 'mapDescriptivesF', 'mapCrossCorrF')
end

toc


%% DX analysis summary at the ML level

maxLayer = 2

MLsummary_XcorrCurvesVelAcf(ML,1,0,'Actin','Vel', maxLayer, ...
    'mapDescriptivesF', 'mapCrossCorrF', 'lagMax0', 40, ...
    'outDirName', 'Xcf_ch1ch0_F_lagMax40')


%% end: DX 



%%  Fluctuation Profiling Around Edge Motion Events (FPEAME)


%%  ML_phaseDescriptives_LB
tic
maxLayer = 2
for i=1:numel(ML.movieDataFile_)
    % Check the tuning parameters of smParam, minimumRunLength, samplingBw within ML_phaseDescriptives function.
    ML_phaseDescriptives(ML, i, maxLayer, 'Actin', 1, 'phaseDescriptives0p4mRL5')

end
toc


%% FPAME summary at the ML level

maxLayer = 2

MLsummary_FluctuationProfiling(ML, 'Actin', maxLayer, 'phaseDescriptives0p4mRL5', ...
    'phaseDesc3_0p4mRL5_lagMax10', 'lagMax0', 10)


%%  end: FPAME



%% ML_meanCV_alongDepth

tic
for i=1:numel(ML.movieDataFile_)
    ML_meanCV_alongDepth(ML, i, 7, 'Actin', 1, 'mapDescriptives_alongDepth')
end
toc


%%  MLsummary_alongDepth

MLsummary_alongDepth(ML,1,7,'Actin','mapDescriptives_alongDepth', ...
     'outDirName', 'topolayersMeanSD3_Chan1')
 
 

%% 
%% for i=1:100; close(figure(i)); end




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
     'omittedWindows', [], 'subFrames', subFr, 'numPerm', 100, 'topograph', 'on', 'parpoolNum', 6)



%%  Channel Descriptives

%iChan = 1;
chanName = chNametag;
chanTitle = chanName;
 
mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, 'adf', 1, ...
     'omittedWindows', [], 'subFrames', subFr, 'numPerm', 100, 'topograph', 'on', 'parpoolNum', 6, 'WithN', 0)

for k=1:300; close(figure(k)); end           

            
%%  specify outputDir (cross correlation analysis)

%figuresDir = fullfile(MD.outputDirectory_, 'mapCrossCorrF');
figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);


%% xcorr chan vs velocity

layerMax = maxLayer;
%
mapXcorrCurvePermutation_Vel(MD, iChan, chanName, layerMax, figuresDir, 'impute', 0, ...
    'subFrames', subFr, 'parpoolNum', 6, 'WithN', 0, 'numPerm', 100, 'lagMax', 40)


for k=1:300; close(figure(k)); end

end


%%
function ML_phaseDescriptives(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
% ML_phaseDescriptives RUN phaseDescriptives_OneChan() at the ML level.
%
% Jungsik Noh, 2018/01/30.

%% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      %ML.getMovies();
pause(1)    


%%  specify outputDir

%figuresDir = fullfile(MD.outputDirectory_, 'phaseDescriptives0p2mRL5_subFr')
figuresDir = fullfile(MD.outputDirectory_, outDirName);


%%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.

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
function ML_meanCV_alongDepth(ML, MDindex, maxLayer, chNametag, iChan, outDirName)
% ML_meanCV_alongDepth RUN mapDescriptives_meanCV_alongDepth() at ML
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


%  omitWin 

%%  

chanName = chNametag;

mapDescriptives_meanCV_alongDepth(MD, iChan, maxLayer, chanName, figuresDir, ...
    'impute',0,'WithN',0, 'subFrames', subFr, 'omittedWindows', [])

for k=1:50; close(figure(k)); end

end


 

%%  EOF
