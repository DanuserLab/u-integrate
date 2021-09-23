function ML_GC_SPAR1ch_informationFlow(ML, MDindex, layerMax, chNametag1, ...
    iChan1, twlagMax, outDirName, varargin)
% ML_iGC_SPAR1ch() RUN MD_iGC_SPAR1ch() for each MD.
%
% Updated:
% J Noh, 2021/02/04. Rename and add comments. 

% get movieData at the MDindex-th location
load(fullfile(ML.movieDataFile_{MDindex}))      %ML.getMovies();
pause(1)                                                %MD = ML.getMovie(MDindex);

disp('===========================')
disp(['======= MDindex = ', num2str(MDindex), ' ======='])
disp('===========================')

load(fullfile(ML.outputDirectory_, 'GCparam.mat'))

ip = inputParser;
ip.addParameter('LB', false);
ip.addParameter('CommonFactorNormAddChVec', {NaN, NaN, NaN});
ip.addParameter('baseOfRatioVec', [NaN, NaN, NaN]);
ip.addParameter('LBoutDirName', 'LBout')
ip.parse(varargin{:});
p = ip.Results;

%%  outputDir

figuresDir = fullfile(MD.outputDirectory_, outDirName)

%%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.

fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
if ~isempty(fInfo)
    load(fullfile(MD.outputDirectory_, 'subFr.mat'));
    disp(['== length of subFr:', num2str(numel(subFr))])
else
    subFr = [];
end

disp(subFr)

%%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
%%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
omitWin = []

if p.LB
    velAnalName = p.LBoutDirName;
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);
    
    omitWin = find(raw.indActive == 0);
    disp(omitWin)
end


%% gc ch1, ch2

chan1Name = chNametag1

MD_GC_SPAR1ch_informationFlow(MD, iChan1, chan1Name, layerMax, ...
    figuresDir, twlagMax, ...
    'WithN', 0, 'omittedWindows', omitWin, 'subFrames', subFr, 'parpoolNum', 10, ...
    'movMedFrameSize', GCparam.movMedFrameSize, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
    'CommonFactorNormAddChVec', p.CommonFactorNormAddChVec, 'factoranMethod', GCparam.factoranMethod, ...
    'baseOfRatioVec', p.baseOfRatioVec, 'EWMA', GCparam.EWMA, ...
    'infoCriterion', GCparam.infoCriterion)

close all

end
