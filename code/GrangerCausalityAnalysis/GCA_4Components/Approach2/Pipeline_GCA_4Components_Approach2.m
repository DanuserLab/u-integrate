%% Main pipeline of Granger-Causality (GC) analysis for 3-channel live cell videos 
%% (Four component analysis pipeline in Approach 2)
%
% UTSW, Dallas, TX (11/2025)
% Authors: Anteneh Godana and Jungsik Noh


%% Workflow
% (Part 1)  
%       Extract subcellular edge motion profiles and molecular activities 
%       from live cell movies, using 'u-register Package'.
% (Part 2 - this pipeline)  
%       Run "Pipeline_GCA_4Components_Approach2.m" to implement Cross correlation analysis, 
%       Fluctuation profiling, and/or Granger-causality analysis. 


%% This main pipeline consists of 6 steps:
%   (1) Set up parameters for Granger-causality analysis (GCA)
%   (2) Quiescent window detection
%   (3) Signal-to-Noise-Ratio (SNR) across different smoothing parameters
%   (4) Activity map visualization, XCorrelation analysis and ADF testing
%   (5) Fluctuation profiling around motion events (FPAME)
%   (6) GCA of 4-components (Approach 2)


%% The only input for this pipeline is the movieList object loaded to workspace.
% Each movieData object in the movieList is required to be already
% processed by segmentation, protrusion profiling, and windowing processes.


%% Step 1. Set up parameters for Granger-causality analysis

GCparam = struct();
% half of the cycle of protrusion/retraction (in frames).
% This needs to be calculated from auto-correlation functions (ACF) of edge
% velocity, which can be calculated by Step 2 below.
% If there is no preliminary information, run Step 2 first to determine this 'halfOfPRCycle' parameter.
halfOfPRCycle = 10;         

GCparam.name1 = 'Actin';    % label of Channel 1 (as listed in the movieData object)
GCparam.name2 = 'Arp3cr';   % label of Channel 2
GCparam.name3 = 'VASP';     % label of Channel 3
GCparam.maxLayer = 4;       % the num of layers to be analyzed

% 2 * half cycle of protrusion/retraction is fed into movMedian period length
% for low frequency normalization.
GCparam.LFNfr = 2*halfOfPRCycle;
% this prefix defines the output folder name
outputName = ['GCA_4comp_LFS_', num2str(GCparam.LFNfr), 'fr_', GCparam.name1, '_', GCparam.name2, '_', GCparam.name3];
GCparam.lagMax = 2*halfOfPRCycle;       % for XCorr analysis
GCparam.movMedFrameSize = 2*halfOfPRCycle;
GCparam.smParamTh = 0.5;
GCparam.factoranMethod = 0000;          
GCparam.infoCriterion = 'AIC';
GCparam.tLagfr = round(halfOfPRCycle/2);    % maximum regression lag in GCA
GCparam.movingAvgSmoothing = false;         % for quiescent window detection
GCparam.EWMA = 0.5;         % lambda=1 (no smoothing), lambda=0.5 (default)

GCparam.outputName = outputName;
disp('==== Defined GC analysis parameters:')
disp(GCparam)
% save the parameter for record in the ML folder
save(fullfile(ML.outputDirectory_, 'GCparam.mat'), 'GCparam')

% define labels of components
Ch1name = GCparam.name1;
LFNCh1name= ['LFS.', GCparam.name1];      % low frequency normalized component 1
LFofCh1name = ['LFof', GCparam.name1];    % low frequency part of component 1
Ch2name = GCparam.name2;
LFNCh2name = ['LFS.', GCparam.name2];
LFofCh2name = ['LFof', GCparam.name2];
Ch3name = GCparam.name3;
LFNCh3name = ['LFS.', GCparam.name3];
LFofCh3name = ['LFof', GCparam.name3];

% define output folder names and key parameters
LBout = 'quiescentWindow_LB';
LFNfr = GCparam.LFNfr;
maxLayer = GCparam.maxLayer;
MDDescOut = ['MapDesc_', outputName];
MDCorrOut = ['MapCCorr_', outputName];

% Channel numbering pattern
%   - 1: channel 1
%   - 21: low frequency subtracted (LFS) channel 1


%% Step 2. Quiescent window detection
%  Identifying quiescent windows by using Ljung-Box test and compute ACF of edge velocities. 

for i=1:numel(ML.movieDataFile_)
    ML_quiescentWindow_Vel_LB(ML, i, LBout, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
        'EWMA', GCparam.EWMA)
end
MLsummary_quiescentWindow(ML, LBout, 'outDirName', ['MLdiagnosticPlots_', outputName])


%% Step 3. Signal-to-Noise-Ratio (SNR) across different smoothing parameters

% ML_SNRoverSmoothParam
for i = 1:numel(ML.movieDataFile_)
    ML_SNRoverSmoothingParams(ML, i, 0, 1, 'Vel', GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')
    ML_SNRoverSmoothingParams(ML, i, 21, maxLayer, LFNCh1name, GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')
    ML_SNRoverSmoothingParams(ML, i, 22, maxLayer, LFNCh2name, GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')
    ML_SNRoverSmoothingParams(ML, i, 23, maxLayer, LFNCh3name, GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')  % New channel
end

% Summary for movie list
MLsummary_SNRoverSmoothingParams(ML, 0, 'Vel', 'SNRoverSmoothParam_EWMA')
MLsummary_SNRoverSmoothingParams(ML, 21, LFNCh1name, 'SNRoverSmoothParam_EWMA')
MLsummary_SNRoverSmoothingParams(ML, 22, LFNCh2name, 'SNRoverSmoothParam_EWMA')
MLsummary_SNRoverSmoothingParams(ML, 23, LFNCh3name, 'SNRoverSmoothParam_EWMA')


%% Step 4. Activity map visualization, XCorrelation analysis and ADF testing

for i=1:numel(ML.movieDataFile_)
    ML_CrossCorr_1chan(ML, i, maxLayer, LFNCh1name, 21, MDDescOut, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_1chan(ML, i, maxLayer, LFNCh2name, 22, MDDescOut, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_1chan(ML, i, maxLayer, LFNCh3name, 23, MDDescOut, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_2chan(ML, i, maxLayer, LFNCh2name, LFNCh1name, 22, 21, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_2chan(ML, i, maxLayer, LFNCh3name, LFNCh1name, 23, 21, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_2chan(ML, i, maxLayer, LFNCh3name, LFNCh2name, 23, 22, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
end

% Summary for the ML
MLsummary_XcorrCurvesVelAcf(ML, 21, 0, LFNCh1name, 'Vel', maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['Xcf_ch21ch0_', outputName])
MLsummary_XcorrCurvesVelAcf(ML, 22, 0, LFNCh2name, 'Vel', maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['Xcf_ch22ch0_', outputName])
MLsummary_XcorrCurvesVelAcf(ML, 23, 0, LFNCh3name, 'Vel', maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['Xcf_ch23ch0_', outputName])
MLsummary_XcorrCurvesVelAcf(ML, 22, 21, LFNCh2name, LFNCh1name, maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['Xcf_ch22ch21_', outputName])
MLsummary_XcorrCurvesVelAcf(ML, 23, 21,LFNCh3name, LFNCh1name,  maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['Xcf_ch23ch21_', outputName])
MLsummary_XcorrCurvesVelAcf(ML, 23, 22, LFNCh3name,LFNCh2name,  maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['Xcf_ch23ch22_', outputName])

% ADF testing
MLsummary_ADFtest(ML, MDDescOut, 0, 'Velocity', 'outDirName', ['MLdiagnosticPlots_', outputName]) 


%% Step 5. Fluctuation profiling around motion events (FPAME)

MDFPOut = ['FluctuationCurves_', outputName, '_Ch1'];
for i = 1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_FPAME(ML, i, maxLayer, LFNCh1name, 21, MDFPOut, 'LB', true, 'LBoutDirName', LBout);
end
% FPAME summary at the ML level
MLsummary_FluctuationProfiling(ML, LFNCh1name, maxLayer, MDFPOut, ...
    ['FluctProfile_', LFNCh1name, '_', outputName], 'lagMax0', GCparam.LFNfr)

MDFPOut2 = ['FluctuationCurves_', outputName, '_Ch2'];
for i = 1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_FPAME(ML, i, maxLayer, LFNCh2name, 22, MDFPOut2, 'LB', true, 'LBoutDirName', LBout);
end
% FPAME summary at the ML level
MLsummary_FluctuationProfiling(ML, LFNCh2name, maxLayer, MDFPOut2, ...
    ['FluctProfile_', LFNCh2name, '_', outputName], 'lagMax0', GCparam.LFNfr)

MDFPOut3 = ['FluctuationCurves_', outputName, '_Ch3'];
for i = 1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_FPAME(ML, i, maxLayer, LFNCh3name, 23, MDFPOut3, 'LB', true, 'LBoutDirName', LBout);
end
% FPAME summary at the ML level
MLsummary_FluctuationProfiling(ML, LFNCh3name, maxLayer, MDFPOut3, ...
    ['FluctProfile_', LFNCh3name, '_', outputName], 'lagMax0', GCparam.LFNfr);


%% Step 6. GCA of 4-components (Approach 2)

MDGCOut = ['GCA_4Components_Approach2_', outputName];
twlagMax0 = [GCparam.tLagfr, 1];

for i=1:numel(ML.movieDataFile_) 
    
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh2name, LFNCh1name,'Vel',LFNCh3name, 22, 21,0,23, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh1name, LFNCh2name,'Vel',LFNCh3name, 21, 22,0,23, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh3name, LFNCh1name,'Vel',LFNCh2name, 23, 21,0,22, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh1name, LFNCh3name,'Vel',LFNCh2name, 21, 23,0,22,...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh2name, LFNCh3name,'Vel',LFNCh1name, 22, 23,0,21, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh3name, LFNCh2name,'Vel',LFNCh1name, 23, 22,0,21, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, 'Vel', LFNCh1name,LFNCh2name,LFNCh3name, 0, 21,22,23, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh1name, 'Vel',LFNCh2name,LFNCh3name, 21, 0,22,23, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, 'Vel', LFNCh2name,LFNCh1name,LFNCh3name, 0, 22,21,23, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh2name, 'Vel',LFNCh1name,LFNCh3name, 22, 0,21,23, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, 'Vel', LFNCh3name,LFNCh1name,LFNCh2name, 0, 23,21,22, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR4ch_Approach2(ML, i, maxLayer, LFNCh3name, 'Vel',LFNCh1name,LFNCh2name, 23, 0,21,22, ...
        twlagMax0, twlagMax0, MDGCOut, 'LB', true, 'LBoutDirName', LBout)
end 

% Summary for ML
MLsummary_iGC_SPAR4ch(ML,LFNCh2name,LFNCh1name,'Vel',LFNCh3name,maxLayer,MDGCOut,'outDirName',MDGCOut)
MLsummary_iGC_SPAR4ch(ML,LFNCh1name,LFNCh2name,'Vel',LFNCh3name,maxLayer,MDGCOut,'outDirName',MDGCOut)

MLsummary_iGC_SPAR4ch(ML,LFNCh3name,LFNCh1name,'Vel',LFNCh2name,maxLayer,MDGCOut,'outDirName',MDGCOut)
MLsummary_iGC_SPAR4ch(ML,LFNCh1name,LFNCh3name,'Vel',LFNCh2name,maxLayer,MDGCOut,'outDirName',MDGCOut)

MLsummary_iGC_SPAR4ch(ML,LFNCh2name,LFNCh3name,'Vel',LFNCh1name,maxLayer,MDGCOut,'outDirName',MDGCOut)
MLsummary_iGC_SPAR4ch(ML,LFNCh3name,LFNCh2name,'Vel',LFNCh1name,maxLayer,MDGCOut,'outDirName',MDGCOut)

MLsummary_iGC_SPAR4ch(ML,'Vel',LFNCh1name,LFNCh2name,LFNCh3name,maxLayer,MDGCOut,'outDirName',MDGCOut)
MLsummary_iGC_SPAR4ch(ML,LFNCh1name,'Vel',LFNCh2name,LFNCh3name,maxLayer,MDGCOut,'outDirName',MDGCOut)

MLsummary_iGC_SPAR4ch(ML,'Vel',LFNCh2name,LFNCh1name,LFNCh3name,maxLayer,MDGCOut,'outDirName',MDGCOut)
MLsummary_iGC_SPAR4ch(ML,LFNCh2name,'Vel',LFNCh1name,LFNCh3name,maxLayer,MDGCOut,'outDirName',MDGCOut)

MLsummary_iGC_SPAR4ch(ML,'Vel',LFNCh3name,LFNCh1name,LFNCh2name,maxLayer,MDGCOut,'outDirName',MDGCOut)
MLsummary_iGC_SPAR4ch(ML,LFNCh3name,'Vel',LFNCh1name,LFNCh2name,maxLayer,MDGCOut,'outDirName',MDGCOut)

%  GC network Graph
disp('==== drawing network diagram')
chanDetailedNames = {LFNCh1name, LFNCh2name, LFNCh3name,  'Vel'};
chanNamesforNet = {GCparam.name1, GCparam.name2, GCparam.name3, 'Edge Velocity'};

GCsummaryDirName = MDGCOut;
GCsummaryDir = fullfile(ML.outputDirectory_, GCsummaryDirName);
maxLayerforNet = GCparam.maxLayer;
% 
GCsummaryToGCnetGraph_4Components(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayerforNet)


%% End of pipeline

disp("== Pipeline_GCA_4Componets_Approach2() is completed.")

%% EOF