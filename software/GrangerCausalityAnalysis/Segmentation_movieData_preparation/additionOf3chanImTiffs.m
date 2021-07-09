%%  2019/05/13
%%  script for additionOf3chanImTiffs

nCh = 3
chname = cell(nCh, 1);
chname{1} = 'Actin-mNG'
chname{2} = 'Arp3cr-Halo'
chname{3} = 'VASP-SNAP'
%
addedChName = 'Actin_Arp3cr_VASP'

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
    outputDir = fullfile(pwd, cellNames{i}, 'Actin_Arp3cr_VASP')
    %
    I = cell(nCh, 1);
    for k = 1:nCh
        fileReads = dir(imgFolderPaths{i}{k});
        ind = arrayfun(@(x) (x.isdir == 0), fileReads);
        
        fileNames = {fileReads(ind).name};
        % only when filename order is not in time order
        fileNames = natsortfiles(fileNames);
        
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
    
    % Addition of k-channels
    I{nCh+1} = uint16(zeros(nImage, mImage, frmax));
    for fr = 1:frmax
        for k = 1:nCh
            I{nCh+1}(:,:,fr) = I{nCh+1}(:,:,fr) + uint16(I{k}(:,:,fr));
        end
    end
    
    if ~isdir(outputDir); mkdir(outputDir); end
    for fr = 1:frmax
        outfname = fullfile(outputDir, [addedChName, '_', sprintf('%04d', fr), '.tif']);
        imwrite(I{nCh+1}(:,:,fr), outfname);
        fprintf(1, '%g ', fr);
    end   
    
end



fprintf('\n')
disp('== end ==')





%%






%%  EOF