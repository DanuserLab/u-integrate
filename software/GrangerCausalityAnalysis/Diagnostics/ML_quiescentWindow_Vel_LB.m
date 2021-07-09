function ML_quiescentWindow_Vel_LB(ML, MDindex, outDirName, varargin)
% ML_quiescentWindow_Vel_LB() RUN mapDescriptives_Vel_LB.m at the movieList
% level to detect quiescent edge segments. %
%
% Updated:
% J Noh, 2021/02/04. Rename and comment. 
% Jungsik Noh, 2018/01/29.

ip = inputParser;
ip.addParameter('movingAvgSmoothing', false);
ip.addParameter('EWMA', 1);

ip.parse(varargin{:});
p = ip.Results;

%%  get movieData

load(fullfile(ML.movieDataFile_{MDindex}));
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

%%  run mapDescriptives_Vel_LB with options

mapDescriptives_Vel_LB(MD, figuresDir, 'impute',1,'movingAvgSmoothing', p.movingAvgSmoothing, ...
    'Folding', 0, 'EWMA', p.EWMA, 'movmeanNum', 5, ...
    'figFlag','off', 'omittedWindows', [], 'subFrames', subFr, 'topograph', 'off')

close all

end
