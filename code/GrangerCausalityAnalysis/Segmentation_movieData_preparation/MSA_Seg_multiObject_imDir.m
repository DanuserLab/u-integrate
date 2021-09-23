function MSA_Seg_multiObject_imDir(inputImgDir, outputDir, varargin)
% MSA_Seg_multiObject_imDir (multi-scale automatic segmentation) 
% Segment a single cell image by combining segmentation
% results obtained at multiple smoothing scales and by three methods 
% (Minmax, Rosin, Otsu). Since it requires only one tuning parameters 
% (tightness) and ‘tightness’=0.5 works well for many cases, 
% it achieves almost automatic segmentation.
%
% Output:
%       If outputDir specifies only a name, the outputDir is generated
%       under the current directory with the specified name. Otherwise, a
%       specific directory can be specified. In the outputDir,
%       a subdirectory named 'MSASeg_masks' will be generated and
%       contain the masks of the images in the inputImgDir.
%
% Usage:    
%       MSA_Seg_multiObject_imDir(inputImgDir, 'MSAmSeg_ch1', 'tightness', 0.5)  
%
% Option:
%       tightness   - Default is 0.5. If 0, the output is the tightest or
%                   smallest masks. If 1, it is the largest masks.
%       numVotes    - A minimum vote scores for a pixel to be the
%                   foreground (0<= score <=42). Default is -1 (not using). 
%                   Specify/Use one of 'tightness' or 'numVotes' threshold. 
%       imagesOut   - Default is true. If true, images overlayed with
%                   boundaries and vote score images (transformed images)
%                   are saved under outputDir.
%
% Updates:
% 2021/07/06. J Noh. Minor edits on footnotes. 
% 2018/05/01, Jungsik Noh. multiObject allowed and Otsu added. VoteScoreImg
% output.
% 2018/04/20, Jungsik Noh. Modified from MSA_Seg_imDir_201802(). Simplify
% options and edit output texts. 
%
% 2017/05/30, Jungsik Noh


ip = inputParser;
ip.addParameter('tightness', 0.5, @(x) isnumeric(x) && (x==-1 || x >= 0 || x<=1));
ip.addParameter('numVotes', -1);
ip.addParameter('imagesOut', 1);
ip.addParameter('figVisible', 'on');
ip.addParameter('finalRefinementRadius', 1);
ip.addParameter('MinimumSize', 10);
ip.addParameter('ObjectNumber', 1000);
%ip.addParameter('parpoolNum', 1);

ip.parse(varargin{:});
p = ip.Results;

if (p.numVotes > 0); p.tightness = -1; end

tic


%% -------- Parameters ---------- %%

masksOutDir = fullfile(outputDir, 'MSASeg_masks');
if ~isdir(masksOutDir); mkdir(masksOutDir); end
    
pString = 'MSA_mask_';      %Prefix for saving masks to file
 

%% Load images

fileReads = dir(inputImgDir);
ind = arrayfun(@(x) (x.isdir == 0), fileReads);

fileNames = {fileReads(ind).name};
frmax = numel(fileNames);

I = cell(frmax, 1);
imgStack = [];
for fr = 1:frmax
    I{fr} = imread(fullfile(inputImgDir, fileNames{fr}));
    imgStack = cat(3, imgStack, I{fr});
    fprintf(1, '%g ', fr);
end


%% TS of 5 numbers

pixelmat = reshape(imgStack, [], frmax);
pixelmat1 = pixelmat;
pixelmat1(pixelmat1 == 0) = NaN;
%sum(isnan(pixelmat1(:)))

mts = mean(pixelmat1, 1, 'omitnan');
medts = median(pixelmat1, 1, 'omitnan');
q1ts = quantile(pixelmat1, 0.25, 1);
q3ts = quantile(pixelmat1, 0.75, 1);
q99ts = quantile(pixelmat1, 0.99, 1);
q01ts = quantile(pixelmat1, 0.01, 1);

fts = figure; 
plot(mts)
hold on

plot(medts)
plot(q1ts)
plot(q3ts)
plot(q01ts)
plot(q99ts)
hold off

legend('Mean', 'Median', 'Perct25', 'Perct75', 'Perct01', 'Perct99')
title('Time series of 5 summary statistics')

%%
saveas(fts, fullfile(outputDir, 'TS_of_5statistics.png'), 'png')
saveas(fts, fullfile(outputDir, 'TS_of_5statistics.fig'), 'fig')


%% due to parfor

%if p.parpoolNum  > 1
%    if isempty(gcp('nocreate')); parpool('local'); end
%end


%% Multi Scale Segmentation

refinedMask = cell(frmax, 1);
voteScoreImgs = cell(frmax, 1); 
currTightness = p.tightness;
currNumVotes = p.numVotes;

    parfor fr = 1:frmax
        disp('=====')
        disp(['Frame: ', num2str(fr)])    
        im = I{fr};
        [refinedMask{fr}, voteScoreImgs{fr}] = multiscaleSeg_multiObject_im(im, ...
            'tightness', currTightness, 'numVotes', currNumVotes, ...
            'finalRefinementRadius', p.finalRefinementRadius, ...
            'MinimumSize', p.MinimumSize, 'ObjectNumber', p.ObjectNumber);

        %Write the refined mask to file
        imwrite(mat2gray(refinedMask{fr}), fullfile(masksOutDir, [pString, fileNames{fr}]) );
    end



%% imagesOut

if p.imagesOut == 1

    if p.numVotes >= 0
        prefname = ['numVotes_', num2str(p.numVotes)];
    elseif p.tightness >= 0
        prefname = ['tightness_', num2str(p.tightness)];
    else 
        prefname = '_';
    end
       
    dName2 = ['MSASeg_maskedImages_', prefname];
    imOutDir = fullfile(outputDir, dName2);
    if ~isdir(imOutDir); mkdir(imOutDir); end

allint = imgStack(:);
intmin = quantile(allint, 0.002);
intmax = quantile(allint, 0.998);

    ftmp = figure('Visible', p.figVisible);
    for fr = 1:frmax
        figure(ftmp)
        imshow(I{fr}, [intmin, intmax])
        hold on
        bdd = bwboundaries(refinedMask{fr});
        
        for k = 1:numel(bdd)
            bdd1 = bdd{k};
            plot(bdd1(:,2), bdd1(:,1), 'r');
        end
        
        %bdd1 = bdd{1};
        %plot(bdd1(:,2), bdd1(:,1), 'r');
        hold off
        
        h = getframe(gcf);
        imwrite(h.cdata, fullfile(imOutDir, fileNames{fr}), 'tif')
    end

%  voteScoreImg
    imOutDir2 = fullfile(outputDir, 'MSASeg_voteScoreImgs');
    if ~isdir(imOutDir2); mkdir(imOutDir2); end   

    for fr = 1:frmax
        imwrite(voteScoreImgs{fr}, fullfile(imOutDir2, ['voteScores_', fileNames{fr}]) );
    end
    
end

%%
toc
disp('Multi-Scale Automatic Segmentation is done!')
disp('==:) close all')

end
