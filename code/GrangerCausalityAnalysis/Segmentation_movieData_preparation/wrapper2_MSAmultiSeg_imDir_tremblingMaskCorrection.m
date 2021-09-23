function wrapper2_MSAmultiSeg_imDir_tremblingMaskCorrection(MSAmSeg_dirName)
%% example_script.m
% for Tada's movies
%%

%inputImgDir =...
%'/project/bioinformatics/Danuser_lab/harvard/raw/Jungsik/crop-Tada-170208-actinmDia1/Cell1/mNG-actin';

%outputDir = ...
%'/project/bioinformatics/Danuser_lab/harvard/raw/Jungsik/crop-Tada-170208-actinmDia1/Cell1/MSASeg-mNG-actin';

% MSAmSeg_dirName = 'MSAmSeg_Actin_Arp3cr_mDia1';
% wrapper2_MSAmultiSeg_imDir_tremblingMaskCorrection('MSAmSeg_Actin-mNG')

%% param
nCh = 1
%chname = cell(nCh, 1);
%chname{nCh} = 'mNG-Arp3' 

%
celllabFolders = dir(pwd);
ind = find([celllabFolders.isdir]);
ind1 = ind(3:end);

cellNames = {celllabFolders(ind1).name}
nCellNames = numel(cellNames)

imgFolderPaths = cell(nCellNames, 1);
for i = 1:nCellNames
    imgFolderPaths{i} = cell(nCh, 1);    % 2 channel
    for k = 1:nCh
        imgFolderPaths{i}{k} = fullfile(pwd, cellNames{i}, MSAmSeg_dirName, 'MSASeg_masks');
    end
end

%imgInd = [3 7 9 10 11 12]

%%
for i = 1:nCellNames

%inputImgDir = fullfile(pwd, 'MSASeg_masks')
% outputDir = fullfile(pwd, 'tremblingCorrected_masks131R5')
inputImgDir = imgFolderPaths{i}{1};
outputDir = fullfile(pwd, cellNames{i}, MSAmSeg_dirName, 'tremblingCorrected_masks131R5')

tremblingMaskCorrectionS131(inputImgDir, outputDir, 'refinementClosureRadius', 5)

%MSA_Seg_imDir(inputImgDir, outputDir, 'type', 'middle', 'parpoolNum', 1)
%MSA_Seg_imDir_201802(inputImgDir, outputDir, 'tightness', 0.5, 'parpoolNum', 12)


end

fprintf('\n')
disp('== end ==')
%%
end
