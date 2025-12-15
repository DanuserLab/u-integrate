function ML_iGC_SPAR4ch_Approach3(ML, MDindex, layerMax, chNametag1, chNametag2, chNametag3, chNametag4, ...
    iChan1, iChan2, iChan3, iChan4, twlagMax, twlagMaxReg, outDirName, varargin)
% ML_iGC_SPAR4ch_update5() RUN MD_iGC_SPAR4ch_updated5() for each MD.
% Anteneh A. Godana update 6/22/2023 addtional channel 
% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}));      % ML.getMovies();
pause(1);                                                % MD = ML.getMovie(MDindex);

disp('===========================');
disp(['======= MDindex = ', num2str(MDindex), ' =======']);
disp('===========================');

load(fullfile(ML.outputDirectory_, 'GCparam.mat'));
%NumLags = 20;
ip = inputParser;
ip.addParameter('LB', false);
ip.addParameter('CommonFactorNormAddChVec', {NaN, NaN, NaN, NaN});
ip.addParameter('baseOfRatioVec', [NaN, NaN, NaN, NaN]);
ip.addParameter('LBoutDirName', 'LBout');
ip.parse(varargin{:});
p = ip.Results;
%NumLags = ip.Results.NumLags;
%% outputDir

figuresDir = fullfile(MD.outputDirectory_, outDirName);

%% Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.

fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
if ~isempty(fInfo)
    load(fullfile(MD.outputDirectory_, 'subFr.mat'));
    disp(['== length of subFr:', num2str(numel(subFr))]);
else
    subFr = [];
end

disp(subFr);

%% outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
omitWin = [];

if p.LB
    velAnalName = p.LBoutDirName;
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0);
    disp(omitWin);
end

%% gc ch1, ch2, ch3, ch4
chan1Name = chNametag1;
chan2Name = chNametag2;
chan3Name = chNametag3;
chan4Name = chNametag4;

MD_iGC_SPAR4ch_Approach3(MD, iChan1, iChan2, iChan3, iChan4, chan1Name, chan2Name, chan3Name, chan4Name, layerMax,  ...
    figuresDir, twlagMax, twlagMaxReg, 'WithN', 0, 'omittedWindows', omitWin, 'subFrames', subFr, 'parpoolNum', 10, ...
    'movMedFrameSize', GCparam.movMedFrameSize, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
    'CommonFactorNormAddChVec', p.CommonFactorNormAddChVec, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatioVec', p.baseOfRatioVec, 'EWMA', GCparam.EWMA, ...
    'infoCriterion', GCparam.infoCriterion);

close all;

end
