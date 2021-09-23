function ML_CrossCorr_2chan(ML, MDindex, layerMax, chNametag1, chNametag2, ...
    iChan1, iChan2, mapXCorrDirName, varargin)
% ML_CrossCorr_2chan() RUN mapXcorrCurvePermutation() for each MD.
%
% Update:
% J Noh, 2021/02/04. Rename and add comments.

%% get movieData at the MDindex-th location

load(fullfile(ML.movieDataFile_{MDindex}))      %ML.getMovies();
pause(1)                                                %MD = ML.getMovie(MDindex);

disp('===========================')
disp(['======= MDindex = ', num2str(MDindex), ' ======='])
disp('===========================')

load(fullfile(ML.outputDirectory_, 'GCparam.mat'))

ip = inputParser;
ip.addParameter('LB', false);
ip.addParameter('CommonFactorNormAddChVec', {NaN, NaN});
ip.addParameter('baseOfRatioVec', [NaN, NaN]);
ip.addParameter('LBoutDirName', 'LBout')
ip.parse(varargin{:});
p = ip.Results;

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
omitWin = []

if p.LB
    %velAnalName = ['mapDescriptives_Vel_LBmv5_', 'EWMA0p5'];
    velAnalName = p.LBoutDirName;
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0);
    disp(omitWin)
end

%% xcorr ch1, ch2

chan1Name = chNametag1
chan2Name = chNametag2

mapXcorrCurvePermutation(MD, iChan1, iChan2, chan1Name, chan2Name,layerMax,figuresDir, 'topograph', 'off', ...
    'impute', 0, 'parpoolNum', 10, 'WithN', 0, 'numPerm', 10, 'omittedWindows', omitWin, 'lagMax', GCparam.lagMax, ...
    'subFrames', subFr, 'movMedFrameSize', GCparam.movMedFrameSize, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
    'CommonFactorNormAddChVec', p.CommonFactorNormAddChVec, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatioVec', p.baseOfRatioVec, 'EWMA', GCparam.EWMA, ...
    'VNtlagMax', GCparam.tLagfr)

close all

end
