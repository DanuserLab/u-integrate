function ML_CrossCorr_1chan(ML, MDindex, maxLayer, chNametag, iChan, ...
    mapDescriptivesDirName, mapXCorrDirName, varargin)
% ML_CrossCorr_1chan() RUN mapDescriptives_OneChan() for the velocity and a
% specified channel and RUN mapXcorrCurvePermutation_Vel() for each MD.
%
% Updated:
% J Noh, 2021/02/04. 'LBoutDirName' added. Rename and add comments. 
% Jungsik Noh, 2019/04/30.

load(fullfile(ML.outputDirectory_, 'GCparam.mat'))

ip = inputParser;
ip.addParameter('LB', false);
ip.addParameter('CommonFactorNormAddCh', NaN);
ip.addParameter('baseOfRatio', 1);
ip.addParameter('LBoutDirName', 'LBout')
ip.parse(varargin{:});
p = ip.Results;

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
omitWin = []

if p.LB
    %velAnalName = ['mapDescriptives_Vel_LBmv5_', 'EWMA0p5'];
    velAnalName = p.LBoutDirName;
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0);
    disp(omitWin)
end

%%  velocity Descriptives

iChan0 = 0;
chan0Name = 'Vel';
chan0Title = 'Velocity (nm/sec)';

mapDescriptives_OneChan(MD, iChan0, 1, chan0Name, chan0Title, figuresDir, 'adf',1,'impute',0, ...
    'movingAvgSmoothing', GCparam.movingAvgSmoothing, 'omittedWindows', omitWin, 'subFrames', subFr,'numPerm', 10, ...
    'parpoolNum', 10, 'topograph', 'off', 'EWMA', GCparam.EWMA, ...
    'VNtlagMax', GCparam.tLagfr)

%%  Chan Descriptives

chanName = chNametag;
chanTitle = chanName;

mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, 'adf', 1, ...
    'omittedWindows', omitWin, 'subFrames', subFr, 'numPerm', 10, 'topograph', 'off', ...
    'parpoolNum', 10, 'WithN', 0, ...
    'movMedFrameSize', GCparam.movMedFrameSize, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio, 'EWMA', GCparam.EWMA, ...
    'VNtlagMax', GCparam.tLagfr)

close all

%%  specify outputDir of correlation analysis

figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);

%% xcorr chan vs velocity

layerMax = maxLayer;
%
mapXcorrCurvePermutation_Vel(MD, iChan, chanName, layerMax, figuresDir, 'impute', 0, ...
    'omittedWindows', omitWin, 'subFrames', subFr, 'parpoolNum', 10, ...
    'WithN', 0, 'numPerm', 10, 'topograph', 'off', 'lagMax', GCparam.lagMax, ...
    'movMedFrameSize', GCparam.movMedFrameSize, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio, 'EWMA', GCparam.EWMA, ...
    'VNtlagMax', GCparam.tLagfr)

close all

end
