function MLsummary_iGC_SPAR4ch(ML, chan1Name, chan2Name, chan3Name, chan4Name, ...
    maxLayer, analName, varargin)
% MLsummary_iGC_SPAR4ch Collect/Summarize GC subcellular P-values computed by
% ML_iGC_SPAR4ch() into per-cell median P-values. Then it constructs GC
% pathway network and connectivity matrices. 
disp('===============================================================')
disp(['== GC 4ch from ', chan1Name, ' to ', chan2Name, ' given ', chan3Name, chan4Name])
disp('===============================================================')

% load ML, example MD
ML.getMovies();
md1 = ML.getMovie(1);
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])
fsaveName0 = ['fr_', chan1Name, '_to_', chan2Name, '_given_', chan3Name,'_' , chan4Name];
disp(['Suffix for the output files is : ', fsaveName0])

%% input parsing

MDs = ML.getMovies();
num = numel(MDs);

ch1ActmapName = [chan1Name];
ch2ActmapName = [chan2Name];
ch3ActmapName = [chan3Name];
ch4ActmapName = [chan4Name];
 
ip = inputParser;
ip.addParameter('outDirName', ['iGC_SPAR4ch_', fsaveName0]);
parse(ip, varargin{:})
p = ip.Results;
fname_GCresult = ['GC_', fsaveName0, '_winFits_Layers.mat'];
%%
disp('================')
disp('num of MovieData')
disp(num)

%% outDir
outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end
%%  GC winFits
maxLayer1= maxLayer-3;
 maxLayer2= maxLayer-2;
 maxLayer3= maxLayer-1;
 maxLayer4= maxLayer-0;
MDs = ML.getMovies();
num = numel(MDs);
%
cellLabels = cell(num, 1);
FpArr = cell(num, maxLayer);
medFpvalMat = nan(num, maxLayer);
medpartR2Mat = nan(num, maxLayer);
medSNRMat = nan(num, maxLayer);
FpArr1 = cell(num, maxLayer1);
medFpvalMat1 = nan(num, maxLayer1);
medpartR2Mat1 = nan(num, maxLayer1);
medSNRMat1 = nan(num, maxLayer1);
medSNRMat1 = nan(num, maxLayer1);
FpArr2 = cell(num, maxLayer2);
medFpvalMat2 = nan(num, maxLayer2);
medpartR2Mat2 = nan(num, maxLayer2);
medSNRMat2 = nan(num, maxLayer2);
medSNRMat2 = nan(num, maxLayer2);
FpArr3 = cell(num, maxLayer3);
medFpvalMat3 = nan(num, maxLayer3);
medpartR2Mat3 = nan(num, maxLayer3);
medSNRMat3 = nan(num, maxLayer3);
medSNRMat3 = nan(num, maxLayer3);
FpArr4 = cell(num, maxLayer4);
medFpvalMat4 = nan(num, maxLayer4);
medpartR2Mat4 = nan(num, maxLayer4);
medSNRMat4 = nan(num, maxLayer4);
medSNRMat4 = nan(num, maxLayer4);
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    cellName = cellLab0(1:end);
    
    load(fullfile(mdDir, analName, fname_GCresult));
    
    S = load(fullfile(mdDir, analName, 'MD_GC_4wSPAR4ch_inputParser.mat'));
    omittedWin = S.p.omittedWindows;
    wmaxi = size(winFits{1}, 1);
    effWinInd = ones(wmaxi, 1);
    effWinInd(omittedWin) = 0;
    
    cellLabels{i} = cellName;

    for indL = 1:maxLayer
        Fpvec = winFits{indL}(:, 8);
        FpArr{i, indL} = Fpvec;
        medFpvalMat(i, indL) = median(Fpvec, 'omitnan');
        partR2vec = winFits{indL}(:, 11);        
        medpartR2Mat(i, indL) = median(partR2vec, 'omitnan');
        SNRvec = winFits{indL}(:, 15);
        medSNRMat(i, indL) = median(SNRvec, 'omitnan');
        
        effSize = sum(~isnan(Fpvec) & effWinInd);
        
        if effSize < 6
            medFpvalMat(i, indL) = nan;
            medpartR2Mat(i, indL) = nan;
            medSNRMat(i, indL) = nan;
        end
    end
for indL1 = 1:maxLayer1
        Fpvec1 = winFits{indL1}(:, 8);
        FpArr1{i, indL1} = Fpvec1;
        medFpvalMat1(i, indL1) = median(Fpvec1, 'omitnan');
        partR2vec1 = winFits{indL1}(:, 11);
        medpartR2Mat1(i, indL1) = median(partR2vec1, 'omitnan');
        SNRvec1 = winFits{indL1}(:, 15);
        medSNRMat1(i, indL1) = median(SNRvec1, 'omitnan');
        
        effSize1 = sum(~isnan(Fpvec1) & effWinInd);
        
        if effSize1 < 6
            medFpvalMat1(i, indL1) = nan;
            medpartR2Mat1(i, indL1) = nan;
            medSNRMat1(i, indL1) = nan;
        end
    end
    for indL2 = 2:maxLayer2
        Fpvec2 = winFits{indL2}(:, 8);
        FpArr2{i, indL2} = Fpvec2;
        medFpvalMat2(i, indL2) = median(Fpvec2, 'omitnan');
        partR2vec2 = winFits{indL2}(:, 11);
        medpartR2Mat2(i, indL2) = median(partR2vec2, 'omitnan');
        SNRvec2 = winFits{indL2}(:, 15);
        medSNRMat2(i, indL2) = median(SNRvec2, 'omitnan');
        
        effSize2 = sum(~isnan(Fpvec2) & effWinInd);
        
        if effSize2 < 6
            medFpvalMat2(i, indL2) = nan;
            medpartR2Mat2(i, indL2) = nan;
            medSNRMat2(i, indL2) = nan;
        end
    end

    for indL3 = 3:maxLayer3
        Fpvec3 = winFits{indL3}(:, 8);
        FpArr3{i, indL3} = Fpvec3;
        medFpvalMat3(i, indL3) = median(Fpvec3, 'omitnan');
        partR2vec3 = winFits{indL3}(:, 11);
        medpartR2Mat3(i, indL3) = median(partR2vec3, 'omitnan');
        SNRvec3 = winFits{indL3}(:, 15);
        medSNRMat3(i, indL3) = median(SNRvec3, 'omitnan');
        
        effSize3 = sum(~isnan(Fpvec3) & effWinInd);
        
        if effSize3 < 6
            medFpvalMat3(i, indL3) = nan;
            medpartR2Mat3(i, indL3) = nan;
            medSNRMat3(i, indL3) = nan;
        end
    end

    for indL4 = 4:maxLayer4
        Fpvec4 = winFits{indL4}(:, 8);
        FpArr4{i, indL4} = Fpvec4;
        medFpvalMat4(i, indL4) = median(Fpvec4, 'omitnan');
        partR2vec4 = winFits{indL4}(:, 11);
        medpartR2Mat4(i, indL4) = median(partR2vec4, 'omitnan');
        SNRvec4 = winFits{indL4}(:, 15);
        medSNRMat4(i, indL4) = median(SNRvec4, 'omitnan');
        
        effSize4 = sum(~isnan(Fpvec4) & effWinInd);
        
        if effSize4 < 6
            medFpvalMat4(i, indL4) = nan;
            medpartR2Mat4(i, indL4) = nan;
            medSNRMat4(i, indL4) = nan;
        end
    end
    end

log10medFpvalMat = log10(medFpvalMat);
log10medFpvalMat1 = log10(medFpvalMat1);
log10medFpvalMat2 = log10(medFpvalMat2);
log10medFpvalMat3 = log10(medFpvalMat3);
log10medFpvalMat4 = log10(medFpvalMat4);

save(fullfile(outDir, ['GC_Fpval_', fsaveName0, '.mat']), ...
    'cellLabels', 'FpArr', 'ch1ActmapName', ...
    'ch2ActmapName','maxLayer', 'analName', 'medFpvalMat', ...
    'log10medFpvalMat', 'medpartR2Mat')

a = 1:maxLayer;
a1 = string(a');
condNames = strcat('Layer', a1)';

tab = array2table(medFpvalMat, 'VariableNames', condNames);
tab2 = array2table(medpartR2Mat, 'VariableNames', condNames);
tabSNR = array2table(medSNRMat, 'VariableNames', condNames);

tabwn = [cell2table(cellLabels), tab];
tab2wn = [cell2table(cellLabels), tab2];
tabSNRwn = [cell2table(cellLabels), tabSNR];

writetable(tabwn, fullfile(outDir, ['GC_medFpval_Table_', fsaveName0, '.csv']))
writetable(tab2wn, fullfile(outDir, ['GC_medpartialR2_Table_', fsaveName0, '.csv']))
writetable(tabSNRwn, fullfile(outDir, ['GC_medSNR_Table_', fsaveName0, '.csv']))

disp(['GC subcellular median Pvalues saved into .csv, ', fullfile(outDir, ['GC_medFpval_Table_', fsaveName0, '.csv'])])
disp(['GC subcellular median partial R-squared saved into .csv, ', fullfile(outDir, ['GC_medpartialR2_Table_', fsaveName0, '.csv'])])
disp(['GC subcellular median SNR saved into .csv, ', fullfile(outDir, ['GC_medSNR_Table_', fsaveName0, '.csv'])])
%% signrank test p-val, 2019/05/16, nan checked
signRankPvec = nan(1, size(medFpvalMat, 2));
for l = 1:size(medFpvalMat, 2)
    pvals = medFpvalMat(:, l);
    if any(~isnan(pvals))
        signRankPvec(l) = signrank(pvals, 0.05, 'tail', 'left');
    else
        signRankPvec(l) = NaN;
    end
end
signRankPvec1 = round(signRankPvec, 3);
medmedFpvalVec = round(median(medFpvalMat, 1, 'omitnan'), 3);

%% significance indicator

sigIndic = num2cell((signRankPvec1 < 0.05)');

%% myBoxplot

fb = figure;
matOut = -log10(medFpvalMat);
matOut = min(6, matOut, 'includenan');

if (size(matOut, 1) > 1)
    f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames);
end

ylim([0, 6.5])

% significance indicator
hold on
for l = 1:size(matOut, 2)
    if sigIndic{l}
        tt = text(l, 6.2, '*', 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'k');
    end
end

h = refline([0 -log10(0.05)]);
h.Color = 'k';
h.LineStyle = '--';
h1 = refline([0 6]);
ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
logptick = -log10(ptick);
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = logptick;
ax.YTickLabel = ptick;
ylabel('Median P-value Per-Cell')
ax.FontSize = 12;

% jittered plot
hold on

myjet = colormap(jet(size(matOut, 1)));

mattmp = 0.1 * randn(size(matOut));
mat2 = mattmp + [1:size(matOut, 2)];

for k = 1:size(matOut, 1)
    s(k) = scatter(mat2(k, :), matOut(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjet(k, :);
end

title0 = ['GC from ', chan1Name, ' to ', chan2Name, ' given ', chan3Name,',', chan4Name];
title2 = ['Median Pval: ', num2str(medmedFpvalVec)];
title3 = ['RankTest Pval: (', num2str(signRankPvec1), ')'];
title({title0; title2; title3})

legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)

saveas(fb, fullfile(outDir, ['GC_BPFpval2_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['GC_BPFpval2_', fsaveName0, '.fig']), 'fig')
saveas(fb, fullfile(outDir, ['GC_BPFpval2_', fsaveName0, '.pdf']), 'pdf')
figure(fb);
title('');
xlabel('');
ylabel('');
ax = gca;
ax.FontSize = 15;
pause(0.1);
saveas3format(fb, outDir, ['GC_BPFpval2_', fsaveName0, '_3format']);
%saveas3format(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '3format']))
%%  myBoxplot

fb = figure('Visible', 'on');

% low bound 3
matOut = 100 * medpartR2Mat;
ymaxR2 = max(10, max(matOut(:)));  % partR2 max is at least 10%

if (size(matOut, 1) > 1)
    f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames);
end
ylim([0, ymaxR2])

ax = gca;
ylabel('Median partial R-squares (%)')
ax.FontSize = 12;

% jittered plot
hold on

myjet = colormap(jet(size(matOut, 1)));

mattmp = 0.05 * randn(size(matOut));
mat2 = mattmp + [1:size(matOut, 2)];

for k = 1:size(matOut, 1)
    s(k) = scatter(mat2(k, :), matOut(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjet(k, :);
end

meanPartR2 = mean(matOut, 1, 'omitnan');
scatter(1:size(matOut, 2), meanPartR2, 120, 'r', '+')

medmedpartR2Mat = round(100 * median(medpartR2Mat, 1, 'omitnan'), 1);
title1 = sprintf("%02.1f%% \t", medmedpartR2Mat(:));

title0 = ['GC partial R-squares from ', chan1Name, ' to ', chan2Name, ' given ', chan3Name,',', chan4Name ];
title({title0; title1})

legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)

saveas(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '.fig']), 'fig')
saveas(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '.pdf']), 'pdf')
figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)
saveas3format(fb, outDir, ['GC_BP_partialRsquare_', fsaveName0, '3format'])

%%
%% Additional plot for SNR

%fc = figure;
fb = figure('Visible', 'on');
matOutSNR = medSNRMat;
%matOutSNR = min(6, matOutSNR, 'includenan');

if (size(matOutSNR, 1) > 1)
    f2 = boxplot(matOutSNR, 'Whisker', Inf, 'Labels', condNames);
end

ylim([min(matOutSNR, [], 'all') - 1, max(matOutSNR, [], 'all') + 1])

hold on

h2 = refline([0 median(matOutSNR, 'all', 'omitnan')]);
h2.Color = 'k';
h2.LineStyle = '--';

ax2 = gca;
ax2.YTickMode = 'manual';
ylabel('Median SNR Per-Cell')
ax2.FontSize = 12;

% jittered plot
hold on

myjetSNR = colormap(jet(size(matOutSNR, 1)));

mattmpSNR = 0.1 * randn(size(matOutSNR));
mat2SNR = mattmpSNR + [1:size(matOutSNR, 2)];
for k = 1:size(matOutSNR, 1)
    s(k) = scatter(mat2SNR(k, :), matOutSNR(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjetSNR(k, :);
end

title0SNR = ['MSNR from ', chan1Name, ' to ', chan2Name, ' given ', chan3Name,', ', chan4Name];
title(title0SNR, 'FontSize', 15)
legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)
saveas3format(fb, outDir, ['GC_percell_boxplot_medianSNR_', fsaveName0, '3format'])
saveas(fb, fullfile(outDir, ['percell_boxplot_medianSNR_', fsaveName0, '.fig']))
saveas(fb, fullfile(outDir, ['percell_boxplot_medianSNR_', fsaveName0, '.png']))
figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)

saveas3format(fb, outDir, ['percell_boxplot_medianSNR_', fsaveName0, '3format'])
%%
%% Additional plot for Median SNR vs Median P-values per cell

fb = figure('Visible', 'on');
log10medFpvalMat = log10(medFpvalMat);
matOut = -log10(medFpvalMat);
matOut = min(6, matOut, 'includenan');
matOutSNR = -log10(medSNRMat);
matOutSNR = min(6, matOutSNR, 'includenan');

h = refline([0 -log10(0.05)]);
h.Color = 'k';
h.LineStyle = '--';
h1 = refline([0 6]);
ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
logptick = -log10(ptick);
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = logptick;
ax.YTickLabel = ptick;
ax.FontSize = 12;
hold on
myjetSNR = colormap(jet(size(matOutSNR, 1)));

for k = 1:size(matOut, 1)
    s(k) = scatter(matOutSNR(k, :), matOut(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjetSNR(k, :);
end
% Axis labels and title
xlabel('Median SNR')
ylabel('Median P-value Per-Cell')

title0SNR = ['MSNR Vs.PVal ', chan1Name, ' to ', chan2Name, ' given ', chan3Name, ', ', chan4Name];
title(title0SNR, 'FontSize', 12)

% Add legend
legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)
% Save plot
saveas(fb, fullfile(outDir, ['SNR_vs_Pval_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['SNR_vs_Pval_', fsaveName0, '.fig']), 'fig')
figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)
saveas3format(fb, outDir, ['SNR_vs_Pval_', fsaveName0, '_3format'])
%% for plot of median P-values vs MSNR for first layer
fb = figure('Visible', 'on');
log10medFpvalMat1 = log10(medFpvalMat1);
matOut1 = -log10(medFpvalMat1);
matOut1 = min(6, matOut1, 'includenan');
matOutSNR1 = -log10(medSNRMat1);
matOutSNR1 = min(6, matOutSNR1, 'includenan');

h = refline([0 -log10(0.05)]);
h.Color = 'k';
h.LineStyle = '--';
h1 = refline([0 6]);
ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
logptick = -log10(ptick);
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = logptick;
ax.YTickLabel = ptick;
ax.FontSize = 12;
% Plot median SNR vs. median P-values
hold on
myjetSNR = colormap(jet(size(matOutSNR1, 1)));

for k = 1:size(matOut, 1)
    s(k) = scatter(matOutSNR1(k, :), matOut1(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjetSNR(k, :);
end
% Axis labels and title
xlabel('Median SNR layer 1')
ylabel('Median P-value Per-Cell layer 1')

title0SNR = ['From ', chan1Name, ' to ', chan2Name, ' given ', chan3Name, ', ', chan4Name];
title(title0SNR, 'FontSize', 15)

% Add legend
legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)
% Save plot
saveas(fb, fullfile(outDir, ['SNR_vs_Pval1_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['SNR_vs_Pval1_', fsaveName0, '.fig']), 'fig')
figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)
saveas3format(fb, outDir, ['SNR_vs_Pval1_', fsaveName0, '_3format'])

%% for plot of median P-values vs MSNR for second layer
fb = figure('Visible', 'on');
log10medFpvalMat2 = log10(medFpvalMat2);
matOut2 = -log10(medFpvalMat2);
matOut2 = min(6, matOut2, 'includenan');
matOutSNR2 = -log10(medSNRMat2);
matOutSNR2 = min(6, matOutSNR2, 'includenan');

h = refline([0 -log10(0.05)]);
h.Color = 'k';
h.LineStyle = '--';
h1 = refline([0 6]);
ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
logptick = -log10(ptick);
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = logptick;
ax.YTickLabel = ptick;
ax.FontSize = 12;
hold on
myjetSNR = colormap(jet(size(matOutSNR2, 1)));

for k = 1:size(matOut, 1)
    s(k) = scatter(matOutSNR2(k, :), matOut2(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjetSNR(k, :);
end
% Axis labels and title
xlabel('Median SNR layer 2')
ylabel('Median P-value Per-Cell layer 2')

title0SNR = ['From', chan1Name, ' to ', chan2Name, ' given ', chan3Name, ', ', chan4Name];
title(title0SNR, 'FontSize', 12)

% Add legend
legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)
% Save plot
saveas(fb, fullfile(outDir, ['SNR_vs_Pval2_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['SNR_vs_Pval2_', fsaveName0, '.fig']), 'fig')
figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)
saveas3format(fb, outDir, ['SNR_vs_Pval2_', fsaveName0, '_3format'])
%% for plot of median P-values vs MSNR for third layer
fb = figure('Visible', 'on');
log10medFpvalMat3 = log10(medFpvalMat3);
matOut3 = -log10(medFpvalMat3);
matOut3 = min(6, matOut3, 'includenan');
matOutSNR3 = -log10(medSNRMat3);
matOutSNR3 = min(6, matOutSNR3, 'includenan');

h = refline([0 -log10(0.05)]);
h.Color = 'k';
h.LineStyle = '--';
h1 = refline([0 6]);
ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
logptick = -log10(ptick);
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = logptick;
ax.YTickLabel = ptick;
%ylabel('Median P-value Per-Cell')
ax.FontSize = 12;
% Plot median SNR vs. median P-values
hold on
myjetSNR = colormap(jet(size(matOutSNR3, 1)));

for k = 1:size(matOut, 1)
    s(k) = scatter(matOutSNR3(k, :), matOut3(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjetSNR(k, :);
end
% Axis labels and title
xlabel('Median SNR layer 3')
ylabel('Median P-value Per-Cell layer 3')

title0SNR = ['From ', chan1Name, ' to ', chan2Name, ' given ', chan3Name, ', ', chan4Name];
title(title0SNR, 'FontSize', 15)

% Add legend
legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)
% Save plot
saveas(fb, fullfile(outDir, ['SNR_vs_Pval3_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['SNR_vs_Pval3_', fsaveName0, '.fig']), 'fig')
figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)
saveas3format(fb, outDir, ['SNR_vs_Pval3_', fsaveName0, '_3format'])

%% for plot of median P-values vs MSNR for fourth layer
fb = figure('Visible', 'on');
log10medFpvalMat4 = log10(medFpvalMat4);
matOut4 = -log10(medFpvalMat4);
matOut4 = min(6, matOut4, 'includenan');
matOutSNR4 = -log10(medSNRMat4);
matOutSNR4 = min(6, matOutSNR4, 'includenan');

h = refline([0 -log10(0.05)]);
h.Color = 'k';
h.LineStyle = '--';
h1 = refline([0 6]);
ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
logptick = -log10(ptick);
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = logptick;
ax.YTickLabel = ptick;
ax.FontSize = 12;
% Plot median SNR vs. median P-values
hold on
myjetSNR = colormap(jet(size(matOutSNR4, 1)));

for k = 1:size(matOut, 1)
    s(k) = scatter(matOutSNR4(k, :), matOut4(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjetSNR(k, :);
end
xlabel('Median SNR layer 4')
ylabel('Median P-value Per-Cell layer 4')

title0SNR = ['From ', chan1Name, ' to ', chan2Name, ' given ', chan3Name, ', ', chan4Name];
title(title0SNR, 'FontSize', 12)

% Add legend
legend(s, cellLabels, 'Location', 'eastoutside', 'FontSize', 5)
% Save plot
saveas(fb, fullfile(outDir, ['SNR_vs_Pval4_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['SNR_vs_Pval4_', fsaveName0, '.fig']), 'fig')
figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)
saveas3format(fb, outDir, ['SNR_vs_Pval4_', fsaveName0, '_3format'])

%% Boxplot for Regression order with Customized Whiskers and Jet Colormap
RegorderMat = nan(num, maxLayer);
 for indL = 1:maxLayer
        allQw0 = winFits{indL}(:, 14);
        FpArr{i, indL} = allQw0 ;
        RegorderMat(i, indL) = median(allQw0,'omitnan') ;    
        effSize = sum(~isnan(allQw0) & effWinInd);
 end
 save(fullfile(outDir, ['Reg_order', fsaveName0, '.mat']), ...
    'cellLabels', 'FpArr', 'ch1ActmapName', ...
    'ch2ActmapName','maxLayer', 'analName', 'medFpvalMat')

tabreg = array2table(RegorderMat, 'VariableNames', condNames);

tabwn = [cell2table(cellLabels), tabreg];

writetable(tabwn, fullfile(outDir, ['Reg_medOrder_Table_', fsaveName0, '.csv']))

disp(['subcellular median reg order saved into .csv,', fullfile(outDir, ['Reg_medOrder_Table_', fsaveName0, '.csv'])])

   matOutro= RegorderMat(i, indL);  
    fb = figure('Visible', 'on');

if (size(matOutro, 1) >1)
    f1 = boxplot(matOutro, 'Whisker', Inf, 'Labels', condNames);
end

hold on

myjet = colormap(jet(size(matOutro, 1)));

mattmp = 0.1 * randn(size(matOutro));
mat2 = mattmp + [1:size(matOutro, 2)];

for k = 1:size(matOutro, 1)
    s(k) = scatter(mat2(k, :), matOutro(k, :), 50);
    s(k).LineWidth = 0.6;
    s(k).MarkerEdgeColor = 'w';
    s(k).MarkerFaceColor = myjet(k, :);
    
end

% Update the title to reflect channel names
title(['Regression order chosen by IC: ' chan1Name ' to ' chan2Name ' given ' chan3Name ', ' chan4Name])
ax = gca;
ax.FontSize = 12;

h.Color = 'k';
h.LineStyle = '--';
h1 = refline([0 6]);
ptick = [0,1, 2, 4, 6, 8 ,10];
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = ptick;
ax.YTickLabel = ptick;

figure(fb);
title(''); xlabel(''); ylabel('')
ax = gca;
ax.FontSize = 15;
pause(0.1)
saveas3format(fb, outDir, ['Regression_order_chosen__', fsaveName0, '3format'])
disp('==== MLsummary_iGC_SPAR4ch is finished!! ====')
end
