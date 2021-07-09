function ML_FPAME(ML, MDindex, maxLayer, chNametag, iChan, outDirName, varargin)
% ML_FPAME RUN phaseDescriptives_OneChan() at the ML level.
%
% Updated:
% J Noh, 2021/02/04. Add 'LBoutDirName'. Rename and add comments. 
% Jungsik Noh, 2018/01/30.

load(fullfile(ML.outputDirectory_, 'GCparam.mat'))

ip = inputParser;
ip.addParameter('LB', false);
ip.addParameter('CommonFactorNormAddCh', NaN);
ip.addParameter('baseOfRatio', NaN);
ip.addParameter('LBoutDirName', 'LBout')
ip.parse(varargin{:});
p = ip.Results;

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
omitWin = []

if p.LB
    %velAnalName = ['mapDescriptives_Vel_LBmv5_', 'EWMA0p5'];
    velAnalName = p.LBoutDirName;
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0)
end

%% phaseMasking to find a proper smoothing parameter

smParamTh = GCparam.smParamTh

[protMask1, retMask1] = phaseMasking(MD, smParamTh, figuresDir, 'minimumRunLength', 3, ...
    'omittedWindows', omitWin, 'subFrames', subFr, ...
    'movingAvgSmoothing', GCparam.movingAvgSmoothing, 'EWMA', GCparam.EWMA);

%%  Retraction / specify sampling bandwidth

Mask0 = retMask1;
%samplingBw = 60
samplingBw = GCparam.lagMax

%%  Chan Descriptives

chanName = [chNametag, '-ret'];
chanTitle = chanName;

%
phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
    'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, ...
    'movingAvgSmoothing', GCparam.movingAvgSmoothing, 'EWMA', GCparam.EWMA, ...
    'movMedFrameSize', GCparam.movMedFrameSize, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio)

%%  Protrusion

Mask0 = protMask1;

%%  Chan Descriptives

chanName = [chNametag, '-prot'];
chanTitle = chanName;

%
phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
    'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, ...
    'movingAvgSmoothing', GCparam.movingAvgSmoothing, 'EWMA', GCparam.EWMA, ...
    'movMedFrameSize', GCparam.movMedFrameSize, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio)

%%  phase MaxMinVel

chanName = chNametag;
chanTitle = chanName;

% minimumRunLength option added
phaseDescriptives_MaxMinVel_OneChan(MD, iChan, maxLayer, chanName, chanTitle, ...
    smParamTh, samplingBw, figuresDir, 'minimumRunLength', 3, ...
    'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, ...
    'movingAvgSmoothing', GCparam.movingAvgSmoothing, 'EWMA', GCparam.EWMA, ...
    'movMedFrameSize', GCparam.movMedFrameSize, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio)

close all

%%
end