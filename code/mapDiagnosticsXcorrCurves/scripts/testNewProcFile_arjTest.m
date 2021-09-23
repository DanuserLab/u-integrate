%% test new process file
%% creat MDs each with 1 channel from raw images 
path1 = '/work/bioinformatics/s170480/Data/Jungsik/mapDDX/To_Develop_DX_FPAME_GUI_forJenny/mDia1-actin-11212013-3-Leedata/';
fullpath1 = fullfile(path1,'images_actin_small');
c1 = Channel(fullpath1);
MD1 = MovieData(c1,path1); % path1 here is MD's outputDirectory_
MD1.setPath(path1); % set movieDataPath_, where to save .mat file

path2 = '/work/bioinformatics/s170480/Data/Jungsik/mapDDX/To_Develop_DX_FPAME_GUI_forJenny/mDia1-actin-12052013-3-LeedataforJenny';
fullpath2 = fullfile(path2,'images_actin_small');
c2 = Channel(fullpath2);
MD2 = MovieData(c2,path2);
MD2.setPath(path2);

MD1.pixelSize_ = 108;
MD1.timeInterval_ = 5;
MD2.pixelSize_ = 108;
MD2.timeInterval_ = 5;

MDs=cell(1,2);
MDs{1} = MD1;
MDs{2} = MD2;
%% setup MDs, add 6 processes
for i = 1:numel(MDs)
    MDs{i}.setFilename(sprintf('moviedata%dactin.mat',i));
    MDs{i}.sanityCheck; % add nFrames_, imSize_, and reader to the object.
    
    % add 6 processes.
    MDs{i}.reset()
    process = MultiScaleAutoSegmentationProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    
    process = MaskRefinementProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    
    process = ProtrusionProcess(MDs{i});
    MDs{i}.addProcess(process);
    funParams = process.funParams_;
    % set 'Mask Refinement' as mask process.
    segInd = MDs{i}.getProcessIndex(MDs{i}.findProcessTag('MaskRefinementProcess'));
    funParams.SegProcessIndex = segInd;
    
    process.setParameters(funParams);
    process.run()
    
    process = WindowingProcess(MDs{i});
    MDs{i}.addProcess(process);
    funParams = process.funParams_;
    % set 'Mask Refinement' as mask process.
    funParams.SegProcessIndex = 2;
    process.setParameters(funParams);
    process.run()
    
    process = ProtrusionSamplingProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    
    process = WindowSamplingProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
end

delete(findall(0)); % close all "java" windows


%% save as a MovieList
ML_dir = '/work/bioinformatics/s170480/Data/Jungsik/mapDDX/To_Develop_DX_FPAME_GUI_forJenny/MLtest_SEP2018/';
ML = MovieList(MDs,ML_dir);
ML.setFilename('movieList.mat')
ML.sanityCheck

ML_path = fullfile(ML_dir,'movieList.mat');
%% save as a MovieList
% ML = MovieList(MDs,'/work/bioinformatics/s184919/Data/Jungsik/Premade_MLMD_2chans/ML2mov');
% ML.setFilename('movieList.mat')
% ML.sanityCheck

%% no need to run above script every time, can load from saved mat files; 
% BUT past run processes all recorded in ML and MD. If need MD with only 6
% previous proceses, then need to re-run the above script.
clear
clc
% load('/work/bioinformatics/s184919/Data/Jungsik/Premade_MLMD_2chans/ML2mov/movieList.mat')
load(ML_path);
% load(ML.movieDataFile_{1})
ML.sanityCheck


MD = ML.movies_{1};
%% precondition / error check (for all 3 processes):
windowingId = MD.getProcessIndex('WindowingProcess');
winProc = MD.getProcess(windowingId);
methodType = winProc.funParams_.MethodName;
if ~isequal(methodType, 'ConstantNumber')
    error("Method used to propagate the windows from one frame to the next was not ConstantNumber")
end
%
winSamplingId = MD.getProcessIndex('WindowSamplingProcess');
if isempty(winSamplingId)
    error("WindowSamplingProcess needs to be done before run this process.")
end
%
if isempty(MD.pixelSize_); error('MD.pixelSize_ is required.'); end
if isempty(MD.timeInterval_); error('MD.timeInterval_ is required.'); end

%%
process = LBTestProcess(MD);
MD.addProcess(process);
process.run()

%%
process = XcorrAnalysisProcess(MD);
MD.addProcess(process);
process.run()

%%
process = FluctuationProfilingProcess(MD);
MD.addProcess(process);
process.run()

%% ML as input testing:
clear 
clc
load('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_0912/ML2mov/movieList.mat')
ML.reset;
% MDs = ML.getMovies();
% MD=MDs{1}
%%
process = LBTestProcessML(ML);
ML.addProcess(process);
process.run()

%
process = XcorrAnalysisProcessML(ML);
ML.addProcess(process);
process.run()
%
process = FluctuationProfilingProcessML(ML);
ML.addProcess(process);
process.run()

%%
package = XcorrFluctuationPackage(ML);
ML.addPackage(package); % addPackage method is in movieObject.m
ML.getPackage(1).createDefaultProcess(1); % getPackage method is in movieObject.m
ML.getPackage(1).getProcess(1).run();
ML.getPackage(1).createDefaultProcess(2); 
ML.getPackage(1).getProcess(2).run();
ML.getPackage(1).createDefaultProcess(3); 
ML.getPackage(1).getProcess(3).run();

%%
% note % {1×1 ThresholdProcess}    {1×1 MaskRefinementProcess} are first 2
% processes of SegmentationPackage.
MD.reset()
package = WindowingPackage(MD); 
MD.addPackage(package);
MD.getPackage(1).createDefaultProcess(1);
MD.getPackage(1).getProcess(1).run();
% MD.getPackage(1).createDefaultProcess(2); 
% MD.getPackage(1).getProcess(2).run();

%% list subfolder
% Get a list of all files and folders in this folder.
files = dir('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_0905/mDia1-actin-11212013-3-Lee data/Xcorr_Analysis')
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)
% Print folder names to command window.
for k = 1 : length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
end

%% list subfolder
d = dir('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_0912/mDia1-actin-11212013-3-Lee data/Xcorr_Analysis');
isub = [d(:).isdir] %# returns logical vector
nameFolds = {d(isub).name}'
nameFolds(ismember(nameFolds,{'.','..'})) = []
