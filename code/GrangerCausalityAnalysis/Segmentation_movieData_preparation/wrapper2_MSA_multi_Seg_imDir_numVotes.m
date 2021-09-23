function wrapper2_MSA_multi_Seg_imDir_numVotes(inputDirName)
%% example_script.m

% from wrapper_MSA_Seg_imDir3()
%% 

%inputImgDir =...
%'/project/bioinformatics/Danuser_lab/harvard/raw/Jungsik/crop-Tada-170208-actinmDia1/Cell1/mNG-actin';

%outputDir = ...
%'/project/bioinformatics/Danuser_lab/harvard/raw/Jungsik/crop-Tada-170208-actinmDia1/Cell1/MSASeg-mNG-actin';
% wrapper2_MSA_multi_Seg_imDir_numVotes('Actin-mNG')

%% param
nCh = 1
chnames = cell(nCh, 1);
%chname{nCh} = 'Actin_Arp3cr_VASP' 
chnames{nCh} = inputDirName 

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
        imgFolderPaths{i}{k} = fullfile(pwd, cellNames{i}, chnames{k});
    end
end

%imgInd = [3 7 9 10 11 12]

%%
for i = 1:nCellNames  

% inputImgDir = fullfile(pwd, 'Actin-mNG')
% outputDir = fullfile(pwd, 'MSAmSeg_Actin-mNG')
inputImgDir = imgFolderPaths{i}{1};
outputDir = fullfile(pwd, cellNames{i}, ['MSAmSeg_numVotes_', chnames{1}])

%MSA_Seg_imDir(inputImgDir, outputDir, 'type', 'middle', 'parpoolNum', 1)
%MSA_Seg_multiObject_imDir(inputImgDir, outputDir, 'tightness', 0.7, ...
%    'ObjectNumber', 1, 'finalRefinementRadius', 5, 'figVisible', 'on')
%MSA_Seg_multiObject_imDir(inputImgDir, outputDir, 'tightness', 0.5, ...
%    'ObjectNumber', 1, 'finalRefinementRadius', 3, 'figVisible', 'on')
%MSA_Seg_multiObject_imDir(inputImgDir, outputDir, 'tightness', 0.7, ...
%    'ObjectNumber', 1, 'finalRefinementRadius', 3, 'figVisible', 'on')
%MSA_Seg_multiObject_imDir(inputImgDir, outputDir, 'numVotes', 22, ...
%    'ObjectNumber', 1, 'finalRefinementRadius', 5, 'figVisible', 'on', 'MinimumSize', 30)
MSA_Seg_multiObject_imDir(inputImgDir, outputDir, 'numVotes', 15, ...
    'ObjectNumber', 1, 'finalRefinementRadius', 5, 'figVisible', 'on')

end


%%
end
