%%  2019/05/13
%%  script for additionOf2chanImTiffs

nCh = 2
chname = cell(nCh, 1);
chname{1} = 'CARMIL1-mNG'
chname{2} = 'Tractin-Halo'

%chname{1} = 'CARMIL1-mNG'
%chname{2} = 'Tractin-Halo'
%chname{1} = 'images_actin'
%chname{2} = 'images_mDia1'
%
addedChName = 'CARMIL1PlusTractin'

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
        imgFolderPaths{i}{k} = fullfile(pwd, cellNames{i}, chname{k});
    end
end

%%  load images
for i = 1:nCellNames
    
    %inputImgDir1 = imgFolderPaths{i}{1};
    %inputImgDir2 = imgFolderPaths{i}{2};
    outputDir = fullfile(pwd, cellNames{i}, addedChName)
    %
    I = cell(nCh, 1);
    for k = 1:nCh
        fileReads = dir(imgFolderPaths{i}{k});
        ind = arrayfun(@(x) (x.isdir == 0), fileReads);
        
        fileNames = {fileReads(ind).name};
        frmax = numel(fileNames);
        
        FileTif= fullfile(imgFolderPaths{i}{k}, fileNames{1});
        InfoImage=imfinfo(FileTif);
        mImage=InfoImage(1).Width
        nImage=InfoImage(1).Height
        
        %I{k} = zeros(nImage, mImage, frmax);
        tmpIm = zeros(nImage, mImage, frmax);
        parfor fr = 1:frmax
            tmpIm(:,:,fr) = imread(fullfile(imgFolderPaths{i}{k}, fileNames{fr}));
            %imgStack = cat(3, imgStack, I{fr});
            %fprintf(1, '%g ', fr);
            disp(fr)
        end
        I{k} = tmpIm;
    end 
    
    % Addition of 2channels
    I{3} = uint16(zeros(nImage, mImage, frmax));
    for fr = 1:frmax
        I{3}(:,:,fr) = I{1}(:,:,fr) + I{2}(:,:,fr);
    end
    
    if ~isdir(outputDir); mkdir(outputDir); end
    for fr = 1:frmax
        outfname = fullfile(outputDir, [addedChName, '_', sprintf('%04d', fr), '.tif']);
        imwrite(I{3}(:,:,fr), outfname);
        fprintf(1, '%g ', fr);
    end   
    
end

disp('==')
disp('== End==')






%%






%%  EOF