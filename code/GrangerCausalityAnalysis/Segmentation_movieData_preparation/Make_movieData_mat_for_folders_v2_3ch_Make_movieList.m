%%  Make_movieData_mat_for_folders
%% 2017/06/23
%%  04/24/2017
%%  02/05/2017

%folderIndex = [3:6, 8:10]

%% Locate in parent folder of image files

% param
nCh = 3
chname = cell(nCh, 1);
chname{1} = 'Actin-mNG'
chname{2} = 'Cytosol-Halo'
chname{3} = 'Ch1-adaptiveRatioOverCh2_bw'

%chname{1} = 'Arp23'
%chname{2} = 'NF2'
%chname{3} = 'pNF2'
%chname{4} = 'Arp23pNF2'   % for seg

%chname{2} = 'EB3disappearance';
%chname{3} = 'GEF_wBlackIm';
%chname{4} = 'GEF_repeated'


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



%%  Locate at the analysis folder

for i =1:nCellNames
%for i = [3 5 7 9 10 11 12]    
    mdDir = fullfile(pwd, cellNames{i})
    if ~isdir(mdDir); mkdir(mdDir); end

    saveFolder = fullfile(mdDir)

    channel = [];
    for k=1:nCh
        tmp = Channel(imgFolderPaths{i}{k});
        channel = [channel, tmp];
    end
    

% Constructor needs an array of channels and an output directory (for analysis)
mdnew = MovieData(channel, saveFolder);

% Set the path where to store the MovieData object.
mdnew.setPath(saveFolder);
mdnew.setFilename('movieData.mat');

% Run sanityCheck on MovieData.
% Check image size and number of frames are consistent.
% Save the movie if successfull
mdnew.sanityCheck; % NOTE I think this sanity check is sometimes annoying

% Set some additional movie properties
%mdnew.numAperture_=1.4;
mdnew.pixelSize_ = 120      % in nm after binning
mdnew.timeInterval_= 3      % in sec


% Save the movieData
mdnew.sanityCheck()
mdnew.save;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Make a movieList
% Start in analysis folder

folderReads = dir(pwd);
ind = find([folderReads.isdir]);
% Param
ind1 = ind(3:end);

folderNames = {folderReads(ind1).name}

%
%folderIndex = [4:11]                            %% input
%folderNames1 = {folderNames{folderIndex}}
fdmax = numel(folderNames)

%  Make a movieList
% Create individual movies

for i = 1 : fdmax
    mdDir = fullfile(pwd, folderNames{i})
    load(fullfile(mdDir, 'movieData.mat'))
    MD
    MDs(i) = MD;
end

% Create movie list
mlDir = fullfile(pwd, 'MLmovies')
if ~isdir(mlDir); mkdir(mlDir); end

ML = MovieList(MDs, mlDir);

% Manipulation via CLI

% Set path properties
ML.setPath(mlDir);
ML.setFilename('movieList.mat');

% Save list
ML.save();
ML.sanityCheck
fprintf(1, 'Movie list saved under: %s\n', ML.getFullPath());


%% Graphical interface

movieSelectorGUI
% Launch viewing interface
%  movieViewer(ML);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EOF
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%