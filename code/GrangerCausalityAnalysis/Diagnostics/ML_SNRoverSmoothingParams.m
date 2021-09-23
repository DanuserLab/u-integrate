function ML_SNRoverSmoothingParams(ML, MDindex, iChan, maxLayer, chanName, arOrder, outDirName, varargin)
% ML_SNRoverSmoothingParams() RUN MD_SNRoverSmoothingParams_EWMA() for each
% movieData of a given movieList. 
%
% Updated:
% J. Noh, 2021/01/27. Rename functions and add comments. 
% Jungsik Noh, 2018/01/29.

%%  get movieData
load(fullfile(ML.movieDataFile_{MDindex}));
pause(1)                  

disp('===========================')
disp(['======= MDindex = ', num2str(MDindex), ' ======='])
disp('===========================')

load(fullfile(ML.outputDirectory_, 'GCparam.mat'))   

ip = inputParser;
ip.addParameter('LB', false);
ip.addParameter('CommonFactorNormAddCh', NaN);
ip.addParameter('baseOfRatio', NaN);
ip.parse(varargin{:});
p = ip.Results;

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

%%  run mapDescriptives_Vel_LB with options

MD_SNRoverSmoothingParams_EWMA(MD, iChan, maxLayer, chanName, arOrder, figuresDir, 'impute',1, 'Folding', 0, ...
    'figFlag','off', 'omittedWindows', [], 'subFrames', subFr, ...
    'movMedFrameSize', GCparam.movMedFrameSize, 'baseOfRatio', p.baseOfRatio)

close all

end
