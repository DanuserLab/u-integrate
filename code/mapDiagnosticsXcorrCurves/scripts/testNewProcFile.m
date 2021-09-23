%% test new process file
%% creat MDs each with 1 channel from raw images 
clear
clc
path1 = '/work/bioinformatics/s184919/Data/Jungsik/copy/mDia1_actin_11212013_3_Lee data';
fullpath1 = fullfile(path1,'images_mDia1');
c1 = Channel(fullpath1);
MD1 = MovieData(c1,path1); % path1 here is MD's outputDirectory_
MD1.setPath(path1); % set movieDataPath_, where to save .mat file

path2 = '/work/bioinformatics/s184919/Data/Jungsik/copy/mDia1_actin_12052013_3_Lee data for Jenny';
fullpath2 = fullfile(path2,'images_mDia1');
c2 = Channel(fullpath2);
MD2 = MovieData(c2,path2);
MD2.setPath(path2);

MDs=cell(1,2);
MDs{1} = MD1;
MDs{2} = MD2;
%% setup MDs, add 6 processes
for i = 1:2
    MDs{i}.setFilename('moviedata%dmDia1.mat');
    MDs{i}.pixelSize_ = 108;
    MDs{i}.timeInterval_ = 5;
    MDs{i}.sanityCheck; % add nFrames_, imSize_, and reader to the object.
    
    % add 6 processes.
    MDs{i}.reset()
    process = MultiScaleAutoSegmentationProcess(MDs{i});
    MDs{i}.addProcess(process);
    process.run()
    
    process = MaskRefinementProcess(MDs{i});
    MDs{i}.addProcess(process);
    funParams = MDs{i}.getProcess(2).funParams_;
    % Manually set Closure Radius parameter
    % Use 5! if use default 3, dataset: mDia1-actin-11212013-3-Lee data/images_mDia1 will not pass protrusion process or windowing process.
    % if use 7, FluctuationProfilingProcess will not pass.
    funParams.ClosureRadius = 5;
    process.setParameters(funParams);
    process.run()
    
    process = ProtrusionProcess(MDs{i});
    MDs{i}.addProcess(process);
    funParams = MDs{i}.getProcess(3).funParams_;
    % set 'Mask Refinement' as mask process.
    funParams.SegProcessIndex = 2;
    process.setParameters(funParams);
    process.run()
    
    process = WindowingProcess(MDs{i});
    MDs{i}.addProcess(process);
    funParams = MDs{i}.getProcess(4).funParams_;
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


% save as a MovieList
ML = MovieList(MDs,'/work/bioinformatics/s184919/Data/Jungsik/copy/ML2mov');
ML.setFilename('movieList1chan.mat')
ML.sanityCheck
%
process = QuiescentWindowDetectionProcessML(ML);
ML.addProcess(process);
process.run()

%
process = ActivityMapDescriptionProcessML(ML);
ML.addProcess(process);
    funParams = ML.getProcess(2).funParams_;
    funParams.omittedWindows = true;   
    funParams.parpoolNum = 6;
    funParams.chanName = {'mDia1'}; 
    process.setParameters(funParams);
process.run()

%% no need to run above script every time, can load from saved mat files; 
clear
clc
load('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_0912/ML2mov/movieList.mat')

load(ML.movieDataFile_{1})
MD

%%%%%%%%%%%
%% creat MDs each with 2 channel from raw images 
clc
clear
path1 = '/work/bioinformatics/s184919/Data/Jungsik/allScript/mDia1_actin_11212013_3_Lee data';
ch1 = fullfile(path1,'images_actin');
chan(1)=Channel(ch1);
ch2 = fullfile(path1,'images_mDia1');
chan(2)=Channel(ch2);
MD1 = MovieData(chan,path1); % path1 here is MD's outputDirectory_
MD1.setPath(path1); % set movieDataPath_, where to save .mat file

path2 = '/work/bioinformatics/s184919/Data/Jungsik/allScript/mDia1_actin_12052013_3_Lee data for Jenny';
ch1 = fullfile(path2,'images_actin');
chans(1)=Channel(ch1);
ch2 = fullfile(path2,'images_mDia1');
chans(2)=Channel(ch2);
MD2 = MovieData(chans,path2); % path1 here is MD's outputDirectory_
MD2.setPath(path2); % set movieDataPath_, where to save .mat file

MDs=cell(1,2);
MDs{1} = MD1;
MDs{2} = MD2;
%% 2nd method to setup MDs, add package and run 6 processes
% made sure using the same parameter set up as GUI.
tic
for i = 1:2
    MDs{i}.setFilename('movieData2chan.mat');
    MDs{i}.pixelSize_ = 108;
    MDs{i}.timeInterval_ = 5;
    MDs{i}.sanityCheck; % add nFrames_, imSize_, and reader to the object.
    
    MDs{i}.reset()
    package = WindowingPackage(MDs{i}); 
    MDs{i}.addPackage(package);
    
    MDs{i}.getPackage(1).createDefaultProcess(1);
    MDs{i}.getPackage(1).getProcess(1).run();
    
    MDs{i}.getPackage(1).createDefaultProcess(2); 
        funParams = MDs{i}.getProcess(2).funParams_;
    % Manually set 'Closure Radius' parameter
    % Use 5! if use default 3, dataset: mDia1-actin-11212013-3-Lee data/images_mDia1 will not pass protrusion process or windowing process.
    % if use 7, FluctuationProfilingProcess will not pass.
        funParams.ClosureRadius = 5;
        MDs{i}.getProcess(2).setParameters(funParams);
    MDs{i}.getPackage(1).getProcess(2).run();
    
    MDs{i}.getPackage(1).createDefaultProcess(3);
        funParams = MDs{i}.getProcess(3).funParams_;
        % set 'Mask Refinement' as mask process.
        funParams.SegProcessIndex = 2;
        funParams.ChannelIndex = 1; % important for this data set.
        MDs{i}.getProcess(3).setParameters(funParams);
    MDs{i}.getPackage(1).getProcess(3).run();
    
    MDs{i}.getPackage(1).createDefaultProcess(4); 
        funParams = MDs{i}.getProcess(4).funParams_;
        % set 'Mask Refinement' as mask process.
        funParams.SegProcessIndex = 2;
        funParams.ChannelIndex = 1; % important for this data set.
        MDs{i}.getProcess(4).setParameters(funParams);
    MDs{i}.getPackage(1).getProcess(4).run();
    
    MDs{i}.getPackage(1).createDefaultProcess(5);
    MDs{i}.getPackage(1).getProcess(5).run();
    
    MDs{i}.getPackage(1).createDefaultProcess(6); 
    MDs{i}.getPackage(1).getProcess(6).run();
    
end

delete(findall(0)); % close all "java" windows
toc % 6.5min on my workstation.
%% save as a MovieList
ML = MovieList(MDs,'/work/bioinformatics/s184919/Data/Jungsik/allScript/ML2mov');
ML.setFilename('movieList2chan.mat')
ML.sanityCheck

%% no need to run above script every time, can load from saved mat files; 
clear
clc
load('/work/bioinformatics/s184919/Data/Jungsik/allScript/ML2mov/movieList2chan.mat')
ML.reset;

%%
tic

package = XcorrFluctuationPackage(ML);
ML.addPackage(package); % addPackage method is in movieObject.m

ML.getPackage(1).createDefaultProcess(1); % getPackage method is in movieObject.m
ML.getPackage(1).getProcess(1).run();

ML.getPackage(1).createDefaultProcess(2); %10/02/2018 changed, default does not use LB Test result.
    % run below 3 lines to use LB Test result:
    funParams = ML.getProcess(2).funParams_;
    funParams.omittedWindows = true;
    funParams.chanName = {'Actin' 'mDia1'}; % rename channel names.
    ML.getProcess(2).setParameters(funParams);   
ML.getPackage(1).getProcess(2).run(); 

ML.getPackage(1).createDefaultProcess(3); %10/02/2018 changed, default does not use LB Test result.
    % run below 3 lines to use LB Test result:
    funParams = ML.getProcess(3).funParams_;
    funParams.omittedWindows = true;
    funParams.chanName = {'Actin' 'mDia1'}; % rename channel names.
    ML.getProcess(3).setParameters(funParams);   
ML.getPackage(1).getProcess(3).run(); 

ML.getPackage(1).createDefaultProcess(4); %10/02/2018 changed, default does not use LB Test result.
    % run below 3 lines to use LB Test result:
    funParams = ML.getProcess(4).funParams_;
    funParams.omittedWindows = true;
    funParams.chanName = {'Actin' 'mDia1'}; % rename channel names.
    ML.getProcess(4).setParameters(funParams);  
ML.getPackage(1).getProcess(4).run();


toc % 21min on BioHPC WebGUI256 RAM; 13 min on my workstation; 17min on WebGPU

%%%%%%%%%%%%%%%

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
process = QuiescentWindowDetectionProcess(MD);
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
load('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_1003/ML2mov/movieList.mat')
ML.reset;
ML1=ML;
load('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_0912/ML2mov/movieList.mat')
ML.reset;
ML2=ML;
% MDs = ML.getMovies();
% MD=MDs{1}

%% test XcorrFluctuationPackage GUI:
% method 1
XcorrFluctuationPackage.GUI(ML)

%% test XcorrFluctuationPackage GUI:
% method 2
XcorrFluctuationPackageGUI(ML)

%%
process = QuiescentWindowDetectionProcessML(ML);
ML.addProcess(process);
process.run()

%%
process = ActivityMapDescriptionProcessML(ML);
ML.addProcess(process);
process.run()
%%
process = XcorrAnalysisProcessML(ML);
ML.addProcess(process); %10/02/2018 changed, default does not use LB Test result.
    % run below 3 lines to use LB Test result:
    funParams = ML.getProcess(2).funParams_;
    funParams.omittedWindows = true;
    process.setParameters(funParams);
process.run()
%%
process = FluctuationProfilingProcessML(ML);
ML.addProcess(process); %10/02/2018 changed, default does not use LB Test result.
    % run below 3 lines to use LB Test result:
    funParams = ML.getProcess(3).funParams_;
    funParams.omittedWindows = true;
    process.setParameters(funParams);
process.run()

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

%% clear subfolders of a dir, no matter it is empty or not
if exist('/work/bioinformatics/s184919/Test2') 
    rmdir('/work/bioinformatics/s184919/Test2','s')
    mkdir('/work/bioinformatics/s184919/Test2')
else 
    mkdir('/work/bioinformatics/s184919/Test2')
end

%% list subfolder
% Get a list of all files and folders in this folder.
files = dir('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_1003/mDia1-actin-12052013-3-Lee data for Jenny/Xcorr_Analysis')
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)
% Print folder names to command window.
for k = 1 : length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
end

%% list subfolder
d = dir('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_1003/mDia1-actin-12052013-3-Lee data/Xcorr_Analysis');
isub = [d(:).isdir] %# returns logical vector
nameFolds = {d(isub).name}'
nameFolds(ismember(nameFolds,{'.','..'})) = []

