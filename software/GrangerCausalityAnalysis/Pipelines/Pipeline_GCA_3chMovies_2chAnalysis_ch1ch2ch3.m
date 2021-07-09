%% /home2/s159113/matlab/GrangerCausalityAnalysis
%% Pipeline_GCA_3chMovies_2chAnalysis_ch1ch2ch3.m

%% The only input for this pipeline is the movieList object loaded to workspace.
% Jungsik Noh, 2018/05/08

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parameters for the pipeline

GCparam = struct();

halfOfPRCycle = 10
GCparam.name1 = 'Actin';
GCparam.name2 = 'Cytosl';
GCparam.name3 = 'ActinNormalized';
GCparam.maxLayer = 4;

% 2*half cycle of protrusion/retraction is fed into movMedian period length
% for low frequency normalization.
GCparam.LFNfr = 2*halfOfPRCycle;
prfxNorm = ['LF', num2str(GCparam.LFNfr), 'fr_PL_GCA_3chMov_2chAnalysis_ewma0p5_', ...
    GCparam.name1, '_', GCparam.name2, '_', GCparam.name3]

GCparam.lagMax = 2*halfOfPRCycle;
GCparam.movMedFrameSize = 2*halfOfPRCycle;
GCparam.smParamTh = 0.5;
GCparam.factoranMethod = 0000;
GCparam.highDimRegMethod = 'IC';
GCparam.infoCriterion = 'AIC';
GCparam.tLagfr = halfOfPRCycle;
GCparam.movingAvgSmoothing = false;
GCparam.EWMA = 0.5;       % lambda = 1 (no smoothing, default)
GCparam.analysisName = 'pipeline_GCA_3chMov-2chAnalysis_ewma0p5';
GCparam.prfxNorm = prfxNorm;

GCparam
save(fullfile(ML.outputDirectory_, 'GCparam.mat'), 'GCparam')
%  end: parameter setup

%% PPL starts with defining chan names
% smParamTh = GCparam.smParamTh

prfxNorm = GCparam.prfxNorm
LFNfr = GCparam.LFNfr
maxLayer = GCparam.maxLayer
 
%
Ezrin = GCparam.name1
LFNEzrin = ['LFS.', GCparam.name1]
LFofEzrin = ['LFof', GCparam.name1]

Arp3 = GCparam.name2
LFNArp3 = ['LFS.', GCparam.name2]
LFofArp3 = ['LFof', GCparam.name2]

mDia1 = GCparam.name3
LFNmDia1 = ['LFS.', GCparam.name3]
LFofmDia1 = ['LFof', GCparam.name3]


%% 2019/04/16, ML_SNRoverSmoothParam

for i=1:numel(ML.movieDataFile_)
    ML_SNRoverSmoothingParams(ML, i, 0, 1, 'Vel', GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')
    ML_SNRoverSmoothingParams(ML, i, 21, maxLayer, LFNEzrin, GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')
    ML_SNRoverSmoothingParams(ML, i, 22, maxLayer, LFNArp3, GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')
    ML_SNRoverSmoothingParams(ML, i, 23, maxLayer, LFNmDia1, GCparam.tLagfr, 'SNRoverSmoothParam_EWMA')
end


MLsummary_SNRoverSmoothingParams(ML, 0, 'Vel', 'SNRoverSmoothParam_EWMA')
MLsummary_SNRoverSmoothingParams(ML, 21, LFNEzrin, 'SNRoverSmoothParam_EWMA')
MLsummary_SNRoverSmoothingParams(ML, 22, LFNArp3, 'SNRoverSmoothParam_EWMA')
MLsummary_SNRoverSmoothingParams(ML, 23, LFNmDia1, 'SNRoverSmoothParam_EWMA')

%%  Diagnotics of vel 
 
%  Identifying quiescent windows by using Ljung-Box test, movmeanNum=5 
LBout = ['quiescentWindow_LB']

for i=1:numel(ML.movieDataFile_)
    ML_quiescentWindow_Vel_LB(ML, i, LBout, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
        'EWMA', GCparam.EWMA)
end
MLsummary_quiescentWindow(ML, LBout, 'outDirName', ['MLdiagnosticPlots_', prfxNorm])

%%  map Diagnostic plots and XCorrelation analysis for all windows

MDDescOut = ['MapDesc_', prfxNorm]   
MDCorrOut = ['MapCCorr_', prfxNorm]

for i=1:numel(ML.movieDataFile_)
% movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName
    ML_CrossCorr_1chan(ML, i, maxLayer, LFNEzrin, 21, MDDescOut, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_1chan(ML, i, maxLayer, LFNArp3, 22, MDDescOut, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_1chan(ML, i, maxLayer, LFNmDia1, 23, MDDescOut, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_2chan(ML, i, maxLayer, LFNArp3, LFNEzrin, 22, 21, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_2chan(ML, i, maxLayer, LFNmDia1, LFNEzrin, 23, 21, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_CrossCorr_2chan(ML, i, maxLayer, LFNmDia1, LFNArp3, 23, 22, MDCorrOut, 'LB', true, 'LBoutDirName', LBout)
end

%% ADF summary after MD descriptives
MLsummary_ADFtest(ML, MDDescOut, 0, 'Velocity', 'outDirName', ['MLdiagnosticPlots_', prfxNorm]) 

MLsummary_XcorrCurvesVelAcf(ML,21,0,LFNEzrin,'Vel', maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, 'outDirName', ['Xcf_ch21ch0_', prfxNorm])
MLsummary_XcorrCurvesVelAcf(ML,22,0,LFNArp3,'Vel', maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, 'outDirName', ['Xcf_ch22ch0_', prfxNorm])
MLsummary_XcorrCurvesVelAcf(ML,23,0,LFNmDia1,'Vel', maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, 'outDirName', ['Xcf_ch23ch0_', prfxNorm])

MLsummary_XcorrCurvesVelAcf(ML,22,21,LFNArp3,LFNEzrin, maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, 'outDirName', ['Xcf_ch22ch21_', prfxNorm])
MLsummary_XcorrCurvesVelAcf(ML,23,21,LFNmDia1,LFNEzrin, maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, 'outDirName', ['Xcf_ch23ch21_', prfxNorm])
MLsummary_XcorrCurvesVelAcf(ML,23,22,LFNmDia1,LFNArp3, maxLayer, ...
    MDDescOut, MDCorrOut, 'lagMax0', LFNfr, 'outDirName', ['Xcf_ch23ch22_', prfxNorm])

%%  ML_phaseDescriptives molecule

MDFPOut = ['FluctuationCurves_', prfxNorm, '_Ch1']
 
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_FPAME(ML, i, maxLayer, LFNEzrin, 21, MDFPOut, 'LB', true, 'LBoutDirName', LBout)
end
%  FPAME summary at the ML level
MLsummary_FluctuationProfiling(ML, LFNEzrin, maxLayer, MDFPOut, ...
    ['FluctProfile_', LFNEzrin, '_0p5mRL3_', prfxNorm], 'lagMax0', GCparam.LFNfr)

%
MDFPOut2 = ['FluctuationCurves_', prfxNorm, '_Ch2']

%maxLayer = 3
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_FPAME(ML, i, maxLayer, LFNArp3, 22, MDFPOut2, 'LB', true, 'LBoutDirName', LBout)
end
%  FPAME summary at the ML level
MLsummary_FluctuationProfiling(ML, LFNArp3, maxLayer, MDFPOut2, ...
    ['FluctProfile_', LFNArp3, '_0p5mRL3_', prfxNorm], 'lagMax0', GCparam.LFNfr)

%
MDFPOut3 = ['FluctuationCurves_', prfxNorm, '_Ch3']

%maxLayer = 3
for i=1:numel(ML.movieDataFile_)
    % check smParam and minimumRunLength within ML_ function.
    ML_FPAME(ML, i, maxLayer, LFNmDia1, 23, MDFPOut3, 'LB', true, 'LBoutDirName', LBout)
end
%  FPAME summary at the ML level
MLsummary_FluctuationProfiling(ML, LFNmDia1, maxLayer, MDFPOut3, ...
    ['FluctProfile_', LFNmDia1, '_0p5mRL3_', prfxNorm], 'lagMax0', GCparam.LFNfr)
 
%%  end: FPAME

%%  GC LFN 
%tic
prfxNorm
MDGCOut2 = ['GCA_2Variables_', prfxNorm, '_lL1wL1tL', num2str(GCparam.tLagfr)]
twlagMax0 = [GCparam.tLagfr, 1]  

for i=1:numel(ML.movieDataFile_) 
    
    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNArp3, LFNEzrin, 22, 21, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNEzrin, LFNArp3, 21, 22, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNEzrin, 'Vel', 21, 0, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR2ch(ML, i, maxLayer, 'Vel', LFNEzrin, 0, 21, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNArp3, 'Vel', 22, 0, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR2ch(ML, i, maxLayer, 'Vel', LFNArp3, 0, 22, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    %
    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNmDia1, LFNEzrin, 23, 21, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNEzrin, LFNmDia1, 21, 23, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNmDia1, 'Vel', 23, 0, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR2ch(ML, i, maxLayer, 'Vel', LFNmDia1, 0, 23, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    %
    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNmDia1, LFNArp3, 23, 22, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR2ch(ML, i, maxLayer, LFNArp3, LFNmDia1, 22, 23, ...
        twlagMax0, twlagMax0, MDGCOut2, 'LB', true, 'LBoutDirName', LBout)   
end
 
 
%% MLsummary GC
%%maxLayer = 3

MLsummary_iGC_SPAR2ch(ML,LFNArp3,LFNEzrin,maxLayer,MDGCOut2,'outDirName',MDGCOut2)
MLsummary_iGC_SPAR2ch(ML,LFNEzrin,LFNArp3,maxLayer,MDGCOut2,'outDirName',MDGCOut2)

MLsummary_iGC_SPAR2ch(ML,LFNEzrin,'Vel',maxLayer,MDGCOut2,'outDirName',MDGCOut2)
MLsummary_iGC_SPAR2ch(ML,'Vel',LFNEzrin,maxLayer,MDGCOut2,'outDirName',MDGCOut2)

MLsummary_iGC_SPAR2ch(ML,LFNArp3,'Vel',maxLayer,MDGCOut2,'outDirName',MDGCOut2)
MLsummary_iGC_SPAR2ch(ML,'Vel',LFNArp3,maxLayer,MDGCOut2,'outDirName',MDGCOut2)

MLsummary_iGC_SPAR2ch(ML,LFNmDia1,LFNEzrin,maxLayer,MDGCOut2,'outDirName',MDGCOut2)
MLsummary_iGC_SPAR2ch(ML,LFNEzrin,LFNmDia1,maxLayer,MDGCOut2,'outDirName',MDGCOut2)

MLsummary_iGC_SPAR2ch(ML,LFNmDia1,'Vel',maxLayer,MDGCOut2,'outDirName',MDGCOut2)
MLsummary_iGC_SPAR2ch(ML,'Vel',LFNmDia1,maxLayer,MDGCOut2,'outDirName',MDGCOut2)

MLsummary_iGC_SPAR2ch(ML,LFNmDia1,LFNArp3,maxLayer,MDGCOut2,'outDirName',MDGCOut2)
MLsummary_iGC_SPAR2ch(ML,LFNArp3,LFNmDia1,maxLayer,MDGCOut2,'outDirName',MDGCOut2)

%%  GCnetwork Graph, GCsummaryToGCnetGraph_merged2ch

disp('==== drawing network diagram')

chan1Name = GCparam.name1;
chan2Name = GCparam.name2;
chanDetailedNames = {LFNEzrin, LFNArp3, 'Vel'};
chanNamesforNet = {GCparam.name1, GCparam.name2, 'Edge Velocity'};

GCsummaryDirName = MDGCOut2
GCsummaryDir = fullfile(ML.outputDirectory_, GCsummaryDirName);
maxLayerforNet = GCparam.maxLayer

GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayerforNet)

%%  GCnetwork Graph, GCsummaryToGCnetGraph_merged2ch

disp('==== drawing network diagram')

chan1Name = GCparam.name1;
chan2Name = GCparam.name3;
chanDetailedNames = {LFNEzrin, LFNmDia1, 'Vel'};
chanNamesforNet = {GCparam.name1, GCparam.name3, 'Edge Velocity'};

GCsummaryDirName = MDGCOut2
GCsummaryDir = fullfile(ML.outputDirectory_, GCsummaryDirName);
maxLayerforNet = GCparam.maxLayer

GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayerforNet)

%%  GCnetwork Graph, GCsummaryToGCnetGraph_merged2ch

disp('==== drawing network diagram')

chan1Name = GCparam.name2;
chan2Name = GCparam.name3;
chanDetailedNames = {LFNArp3, LFNmDia1, 'Vel'};
chanNamesforNet = {GCparam.name2, GCparam.name3, 'Edge Velocity'};

GCsummaryDirName = MDGCOut2
GCsummaryDir = fullfile(ML.outputDirectory_, GCsummaryDirName);
maxLayerforNet = GCparam.maxLayer

GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayerforNet)


%%  GC LFN, ML_GC2_SPAR3ch
%tic
prfxNorm
MDGCOut3 = ['GCA_ch1ch2_3Variables_', prfxNorm, '_lL1wL1tL', num2str(GCparam.tLagfr)]
twlagMax0 = [GCparam.tLagfr, 1] 

for i=1:numel(ML.movieDataFile_) 
    
    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNArp3, LFNEzrin,'Vel', 22, 21,0, ...
        twlagMax0, twlagMax0, MDGCOut3, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNEzrin, LFNArp3,'Vel', 21, 22,0, ...
        twlagMax0, twlagMax0, MDGCOut3, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNEzrin, 'Vel',LFNArp3, 21, 0,22, ...
        twlagMax0, twlagMax0, MDGCOut3, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, 'Vel', LFNEzrin,LFNArp3, 0, 21,22, ...
        twlagMax0, twlagMax0, MDGCOut3, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNArp3, 'Vel',LFNEzrin, 22, 0,21, ...
        twlagMax0, twlagMax0, MDGCOut3, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, 'Vel', LFNArp3,LFNEzrin, 0, 22,21, ...
        twlagMax0, twlagMax0, MDGCOut3, 'LB', true, 'LBoutDirName', LBout)
end 
 
MLsummary_iGC_SPAR3ch(ML,LFNArp3,LFNEzrin,'Vel',maxLayer,MDGCOut3,'outDirName',MDGCOut3)
MLsummary_iGC_SPAR3ch(ML,LFNEzrin,LFNArp3,'Vel',maxLayer,MDGCOut3,'outDirName',MDGCOut3)

MLsummary_iGC_SPAR3ch(ML,LFNEzrin,'Vel',LFNArp3,maxLayer,MDGCOut3,'outDirName',MDGCOut3)
MLsummary_iGC_SPAR3ch(ML,'Vel',LFNEzrin,LFNArp3,maxLayer,MDGCOut3,'outDirName',MDGCOut3)

MLsummary_iGC_SPAR3ch(ML,LFNArp3,'Vel',LFNEzrin,maxLayer,MDGCOut3,'outDirName',MDGCOut3)
MLsummary_iGC_SPAR3ch(ML,'Vel',LFNArp3,LFNEzrin,maxLayer,MDGCOut3,'outDirName',MDGCOut3)

%%  GCnetwork Graph, GCsummaryToGCnetGraph_merged2ch

disp('==== drawing network diagram')

chan1Name = GCparam.name1;
chan2Name = GCparam.name2;
chanDetailedNames = {LFNEzrin, LFNArp3, 'Vel'};
chanNamesforNet = {GCparam.name1, GCparam.name2, 'Edge Velocity'};

GCsummaryDirName = MDGCOut3
GCsummaryDir = fullfile(ML.outputDirectory_, GCsummaryDirName);
maxLayerforNet = GCparam.maxLayer

GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayerforNet)


%%
prfxNorm
MDGCOut4 = ['GCA_ch1ch3_3Variables_', prfxNorm, '_lL1wL1tL', num2str(GCparam.tLagfr)]
twlagMax0 = [GCparam.tLagfr, 1]
%twlagMax0 = [6, 1] 

%maxLayer = 3 
% i=4   
for i=1:numel(ML.movieDataFile_)
    % outputDir = MD_2chan_2factors_Xcorr
    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNmDia1, LFNEzrin,'Vel', 23, 21,0, ...
        twlagMax0, twlagMax0, MDGCOut4, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNEzrin, LFNmDia1,'Vel', 21, 23,0, ...
        twlagMax0, twlagMax0, MDGCOut4, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNEzrin, 'Vel',LFNmDia1, 21, 0,23, ...
        twlagMax0, twlagMax0, MDGCOut4, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, 'Vel', LFNEzrin,LFNmDia1, 0, 21,23, ...
        twlagMax0, twlagMax0, MDGCOut4, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNmDia1, 'Vel',LFNEzrin, 23, 0,21, ...
        twlagMax0, twlagMax0, MDGCOut4, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, 'Vel', LFNmDia1,LFNEzrin, 0, 23,21, ...
        twlagMax0, twlagMax0, MDGCOut4, 'LB', true, 'LBoutDirName', LBout)
end 
 
MLsummary_iGC_SPAR3ch(ML,LFNmDia1,LFNEzrin,'Vel',maxLayer,MDGCOut4,'outDirName',MDGCOut4)
MLsummary_iGC_SPAR3ch(ML,LFNEzrin,LFNmDia1,'Vel',maxLayer,MDGCOut4,'outDirName',MDGCOut4)

MLsummary_iGC_SPAR3ch(ML,LFNEzrin,'Vel',LFNmDia1,maxLayer,MDGCOut4,'outDirName',MDGCOut4)
MLsummary_iGC_SPAR3ch(ML,'Vel',LFNEzrin,LFNmDia1,maxLayer,MDGCOut4,'outDirName',MDGCOut4)

MLsummary_iGC_SPAR3ch(ML,LFNmDia1,'Vel',LFNEzrin,maxLayer,MDGCOut4,'outDirName',MDGCOut4)
MLsummary_iGC_SPAR3ch(ML,'Vel',LFNmDia1,LFNEzrin,maxLayer,MDGCOut4,'outDirName',MDGCOut4)

%%  GCnetwork Graph, GCsummaryToGCnetGraph_merged2ch

disp('==== drawing network diagram')

chan1Name = GCparam.name1;
chan2Name = GCparam.name3;
chanDetailedNames = {LFNEzrin, LFNmDia1, 'Vel'};
chanNamesforNet = {GCparam.name1, GCparam.name3, 'Edge Velocity'};

GCsummaryDirName = MDGCOut4
GCsummaryDir = fullfile(ML.outputDirectory_, GCsummaryDirName);
maxLayerforNet = GCparam.maxLayer

GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayerforNet)


%%
prfxNorm
MDGCOut5 = ['GCA_ch2ch3_3Variables_', prfxNorm, '_lL1wL1tL', num2str(GCparam.tLagfr)]
twlagMax0 = [GCparam.tLagfr, 1]
%twlagMax0 = [6, 1] 

%maxLayer = 3 
% i=4   
for i=1:numel(ML.movieDataFile_)
    % outputDir = MD_2chan_2factors_Xcorr
    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNmDia1, LFNArp3,'Vel', 23, 22,0, ...
        twlagMax0, twlagMax0, MDGCOut5, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNArp3, LFNmDia1,'Vel', 22, 23,0, ...
        twlagMax0, twlagMax0, MDGCOut5, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNArp3, 'Vel',LFNmDia1, 22, 0,23, ...
        twlagMax0, twlagMax0, MDGCOut5, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, 'Vel', LFNArp3,LFNmDia1, 0, 22,23, ...
        twlagMax0, twlagMax0, MDGCOut5, 'LB', true, 'LBoutDirName', LBout)

    ML_iGC_SPAR3ch(ML, i, maxLayer, LFNmDia1, 'Vel',LFNArp3, 23, 0,22, ...
        twlagMax0, twlagMax0, MDGCOut5, 'LB', true, 'LBoutDirName', LBout)
    ML_iGC_SPAR3ch(ML, i, maxLayer, 'Vel', LFNmDia1,LFNArp3, 0, 23,22, ...
        twlagMax0, twlagMax0, MDGCOut5, 'LB', true, 'LBoutDirName', LBout)
end 
 
MLsummary_iGC_SPAR3ch(ML,LFNmDia1,LFNArp3,'Vel',maxLayer,MDGCOut5,'outDirName',MDGCOut5)
MLsummary_iGC_SPAR3ch(ML,LFNArp3,LFNmDia1,'Vel',maxLayer,MDGCOut5,'outDirName',MDGCOut5)

MLsummary_iGC_SPAR3ch(ML,LFNArp3,'Vel',LFNmDia1,maxLayer,MDGCOut5,'outDirName',MDGCOut5)
MLsummary_iGC_SPAR3ch(ML,'Vel',LFNArp3,LFNmDia1,maxLayer,MDGCOut5,'outDirName',MDGCOut5)

MLsummary_iGC_SPAR3ch(ML,LFNmDia1,'Vel',LFNArp3,maxLayer,MDGCOut5,'outDirName',MDGCOut5)
MLsummary_iGC_SPAR3ch(ML,'Vel',LFNmDia1,LFNArp3,maxLayer,MDGCOut5,'outDirName',MDGCOut5)

%%  GCnetwork Graph, GCsummaryToGCnetGraph_merged2ch

disp('==== drawing network diagram')

chan1Name = GCparam.name2;
chan2Name = GCparam.name3;
chanDetailedNames = {LFNArp3, LFNmDia1, 'Vel'};
chanNamesforNet = {GCparam.name2, GCparam.name3, 'Edge Velocity'};

GCsummaryDirName = MDGCOut5
GCsummaryDir = fullfile(ML.outputDirectory_, GCsummaryDirName);
maxLayerforNet = GCparam.maxLayer

GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayerforNet)


%%%%
%%  GC_intra  toward OILR

prfxNorm
MDGCintraOut = ['GCintra', '_', prfxNorm, '_lL1wL1tL', num2str(GCparam.tLagfr)]
twlagMax0 = [GCparam.tLagfr, 1]

for i=1:numel(ML.movieDataFile_) 
    
    ML_GC_SPAR1ch_informationFlow(ML, i, maxLayer, LFNEzrin, 21, ...
        twlagMax0, MDGCintraOut, 'LB', true, 'LBoutDirName', LBout)
    ML_GC_SPAR1ch_informationFlow(ML, i, maxLayer, LFNArp3, 22, ...
        twlagMax0, MDGCintraOut, 'LB', true, 'LBoutDirName', LBout)
    ML_GC_SPAR1ch_informationFlow(ML, i, maxLayer, LFNmDia1, 23, ...
        twlagMax0, MDGCintraOut, 'LB', true, 'LBoutDirName', LBout)
end

%  GCintra summary 
MLsummary_GC_SPAR1ch_informationFlow(ML,21,LFNEzrin,maxLayer,MDGCintraOut,'outDirName',MDGCintraOut)
MLsummary_GC_SPAR1ch_informationFlow(ML,22,LFNArp3,maxLayer,MDGCintraOut,'outDirName',MDGCintraOut)
MLsummary_GC_SPAR1ch_informationFlow(ML,23,LFNmDia1,maxLayer,MDGCintraOut,'outDirName',MDGCintraOut)

%% 2019/12/12 LayerCC  
MDfieldXcorrOut = ['InterLayerCCorr_', prfxNorm]

for i=1:numel(ML.movieDataFile_) 
    
    ML_layerCrossCorr(ML, i, maxLayer, LFNEzrin, LFNEzrin, 21, 21, MDfieldXcorrOut, 'LB', true, 'LBoutDirName', LBout)
    ML_layerCrossCorr(ML, i, maxLayer, LFNArp3, LFNArp3, 22, 22, MDfieldXcorrOut, 'LB', true, 'LBoutDirName', LBout)        
    ML_layerCrossCorr(ML, i, maxLayer, LFNmDia1, LFNmDia1, 23, 23, MDfieldXcorrOut, 'LB', true, 'LBoutDirName', LBout)        
end

%% fieldXcorr summary at the ML level 

MLsummary_layerCrossCorrCurves(ML,21,21,LFNEzrin,LFNEzrin, maxLayer, MDfieldXcorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['InterLayerXcf_ch21_', prfxNorm])
MLsummary_layerCrossCorrCurves(ML,22,22,LFNArp3,LFNArp3, maxLayer, MDfieldXcorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['InterLayerXcf_ch22_', prfxNorm])
MLsummary_layerCrossCorrCurves(ML,23,23,LFNmDia1,LFNmDia1, maxLayer, MDfieldXcorrOut, 'lagMax0', LFNfr, ...
    'outDirName', ['InterLayerXcf_ch23_', prfxNorm])

%
disp('== End of Pipeline_GCA_3chMovies_2chAnalysis_ch1ch2ch3 ==')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2020/01/17, IntraGC Visual Presentation
%prfxNorm
%MDOut_visualEdgeOverTime = ['VisualEdgeOverTime', '_', prfxNorm, '_lL1wL1tL', num2str(GCparam.tLagfr)]





%%  EOF
