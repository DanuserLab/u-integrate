function MLsummary_iGC_SPAR2ch(ML, chan1Name, chan2Name, ...
    maxLayer, analName, varargin)
% MLsummary_iGC_SPAR2ch Collect/Summarize GC subcellular P-values computed by
% ML_iGC_SPAR2ch() into per-cell median P-values. Then it constructs GC
% pathway network and connectivity matrices. 
%
% Option:
%       outDirName  - Specify a name of output directory.
%
% Updates:
% J Noh, 2021/02/05. Rename and add comments. 
% J Noh, 2019/12/27. Add partial Rsquare summary. Remove ttest threshold.
% J Noh, 2017/09/25. Include the summary of ACF of channels.
% Jungsik Noh, 2017/05/17

disp('===============================================================')
disp(['== GC 2ch from ', chan1Name, ' to ', chan2Name])
disp('===============================================================')

% load ML, example MD
ML.getMovies();
md1 = ML.getMovie(1);
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])

%fsaveName0 = ['frCh', num2str(iChan1), 'toCh', num2str(iChan2)];
fsaveName0 = ['fr_', chan1Name, '_to_', chan2Name];

disp(['Suffix for the output files is : ', fsaveName0])

%% input parsing

MDs = ML.getMovies();
num = numel(MDs);

ch1ActmapName = [chan1Name];
ch2ActmapName = [chan2Name];
 
ip = inputParser;
ip.addParameter('outDirName', ['iGC_SPAR2ch_', fsaveName0]);
%ip.addParameter('timeInterval', md1.timeInterval_);
%ip.addParameter('ccfLagMax', -1);
parse(ip, varargin{:})
p = ip.Results;

fname_GCresult = ['GC_', fsaveName0, '_winFits_Layers.mat'];
%MDtimeIntvl = p.timeInterval;

%%
disp('================')
disp('num of MovieData')
disp(num)

%% outDir

outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%%  GC winFits

MDs = ML.getMovies();
num = numel(MDs);
%
cellLabels = cell(num, 1);
FpArr = cell(num, maxLayer);
medFpvalMat = nan(num, maxLayer);
medpartR2Mat = nan(num, maxLayer);

for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    cellName = cellLab0(1:end);    % Generic cellName
    %
    load(fullfile(mdDir, analName, fname_GCresult));    
    % 2019/06
    S = load(fullfile(mdDir, analName, 'MD_GC_4wSPAR3ch_inputParser.mat'));
    omittedWin = S.p.omittedWindows;
    wmaxi = size(winFits{1}, 1);
    effWinInd = ones(wmaxi, 1);
    effWinInd(omittedWin) = 0;
    
    cellLabels{i} = cellName;
    %copyfile(fullfile(mdDir, analName, ['PW_xcCurve_', ch12fname, '.png']), ...
    %    fullfile(outDir, [cellLabels{i}, '_PW_xcorrCurve_', ch12fname, '.png']) )
    
    for indL = 1:maxLayer
        Fpvec = winFits{indL}(:, 8);
        FpArr{i, indL} = Fpvec;
        medFpvalMat(i, indL) = median(Fpvec, 'omitnan');    % median P-values
        partR2vec = winFits{indL}(:, 11);        
        medpartR2Mat(i, indL) = median(partR2vec, 'omitnan'); % median partial R-square
        
        effSize = sum(~isnan(Fpvec) & effWinInd);
        %disp(effSize)
        if effSize <  6                %0.5*sum(effWinInd)
            medFpvalMat(i, indL) = nan;
            medpartR2Mat(i, indL) = nan;
        end
    end
    
end

log10medFpvalMat = log10(medFpvalMat);
%
save(fullfile(outDir, ['GC_Fpval_', fsaveName0, '.mat']), ...
    'cellLabels', 'FpArr', 'ch1ActmapName', ...
    'ch2ActmapName', 'maxLayer', 'analName', 'medFpvalMat', ...
    'log10medFpvalMat', 'medpartR2Mat')

%% writetable medFpvalMat

a = 1:maxLayer;
%a1 = cellstr(num2str(a'));
a1 = string(a');
condNames = strcat('Layer', a1)';

tab = array2table(medFpvalMat, 'VariableNames', condNames);
tabwn = [cell2table(cellLabels), tab];
tab2 = array2table(medpartR2Mat, 'VariableNames', condNames);
tab2wn = [cell2table(cellLabels), tab2];

disp(tabwn)

writetable(tabwn, fullfile(outDir, ['GC_medFpval_Table_', fsaveName0, '.csv']))
writetable(tab2wn, fullfile(outDir, ['GC_medpartialR2_Table_', fsaveName0, '.csv']))

%%  signrank test p-val, 2019/05/16, nan checked

signRankPvec = nan(1, size(medFpvalMat,2));
for l = 1:size(medFpvalMat, 2)
    pvals = medFpvalMat(:, l);
    if any(~isnan(pvals));
        signRankPvec(l) = signrank(pvals, 0.05, 'tail', 'left');
    else
        signRankPvec(l) = NaN;
    end
end
signRankPvec1 = round(signRankPvec, 3);
medmedFpvalVec = round(median(medFpvalMat, 1, 'omitnan'), 3);

%% significance indicator

sigIndic = num2cell((signRankPvec1 < 0.05)');

%%  myBoxplot

fb = figure;
%matOut = max(-3, -log10medFpvalMat);
matOut = -log10medFpvalMat;
matOut = min(6, matOut, 'includenan');
%
if (size(matOut, 1) > 1)
    f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames);
end
%ylim([0, max(matOut(:))*1.1])
ylim([0, 6.5])
% significance indicator
hold on
for l = 1:size(matOut, 2)
   if sigIndic{l}
       tt = text(l, 6.2, '*', 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'k'); 
   end
end

h = refline([0 -log10(0.05)]);h.Color='k'; h.LineStyle='--';
h1 = refline([0 6]);  %h.Color='k'; h.LineStyle='--';
ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
logptick = -log10(ptick);
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = logptick;
ax.YTickLabel = ptick;    %{'0.05', '0.01', '0.001', ' '};
ylabel('Median P-value Per-Cell')

ax.FontSize = 12;
%ax.XTickLabelRotation = 15;

% jittered plot
hold on

myjet = colormap(jet(size(matOut, 1)));

mattmp = 0.1*randn(size(matOut));
mat2 = mattmp + [1:size(matOut, 2)];

for k = 1:size(matOut, 1)
    s(k) = scatter(mat2(k,:), matOut(k,:), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjet(k,:);
    %s.MarkerFaceColor = 'b';
end

title0 = ['GC from ', chan1Name, ' to ', chan2Name ];
%title1 = ['from ', chan1Name, ' to ', chan2Name, ' given ', chan3Name];
title2 = ['Median Pval : ', num2str(medmedFpvalVec)];
title3 = ['RankTest Pval: (', num2str(signRankPvec1), ')'];
title({title0;  title2; title3})

legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)

%%
saveas(fb, fullfile(outDir, ['GC_BPFpval2_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['GC_BPFpval2_', fsaveName0, '.fig']), 'fig')

%%

figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca; ax.FontSize = 15;
%ax.XTickLabelRotation = 15;
pause(0.1)

saveas3format(fb, outDir, ['GC_BPFpval2_', fsaveName0, '_3format'])

%%  myBoxplot, medpartR2Mat

fb = figure('Visible', 'on');

% low bound 3
matOut = 100 * medpartR2Mat;
ymaxR2 = max(10, max(matOut(:)));  % partR2 max is at least 10%

%
if (size(matOut, 1) > 1)
    f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames);
end
%ylim([0, max(matOut(:))*1.1])
ylim([0, ymaxR2])

%h = refline([0 -log10(0.05)]);h.Color='k'; h.LineStyle='--';
%h1 = refline([0 6]);  %h.Color='k'; h.LineStyle='--';
%ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
%logptick = -log10(ptick);
ax = gca;
%ax.YTickMode = 'manual';
%ax.YTick = logptick;
%ax.YTickLabel = ptick;    %{'0.05', '0.01', '0.001', ' '};
ylabel('Median partial R-squares (%)')

ax.FontSize = 12;

% jittered plot
hold on

myjet = colormap(jet(size(matOut, 1)));

mattmp = 0.05*randn(size(matOut));
mat2 = mattmp + [1:size(matOut, 2)];

for k = 1:size(matOut, 1)
    s(k) = scatter(mat2(k,:), matOut(k,:), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjet(k,:);
    %s.MarkerFaceColor = 'b';
end

%s = scatter(mat2(:), matOut(:), [], 1:size(matOut, 1));
%s.MarkerFaceColor = [0 0.5 0.5];
% -log10(averaged P-value)
%meanPval = mean(medFpvalMat, 1, 'omitnan');
meanPartR2 = mean(matOut, 1, 'omitnan');
scatter(1:size(matOut, 2), meanPartR2, 120, 'r', '+')
medmedpartR2Mat = round(100*median(medpartR2Mat, 1, 'omitnan'), 1);
title1 = sprintf("%02.1f%% \t", medmedpartR2Mat(:));

title0 = ['GC partial R-squares from ', chan1Name, ' to ', chan2Name ];
%
%title1 = num2str(medmedpartR2Mat);
title({title0; title1})

legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)


%%
saveas(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '.fig']), 'fig')

%%

figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca; ax.FontSize = 15;
%ax.XTickLabelRotation = 15;
pause(0.1)

saveas3format(fb, outDir, ['GC_BP_partialRsquare_', fsaveName0, '_3format'])

%%
disp('==== MLsummary_iGC_SPAR2ch is finished!! ====')
disp('== :) close all ')

end
