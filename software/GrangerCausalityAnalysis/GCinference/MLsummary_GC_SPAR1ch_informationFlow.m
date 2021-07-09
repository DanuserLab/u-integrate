function MLsummary_GC_SPAR1ch_informationFlow(ML, iChan1, chan1Name, ...
    maxLayer, analName, varargin)
% MLsummary_GC_SPAR1ch_informationFlow Collect/Summarize GC subcellular P-values computed by
% ML_GC_SPAR1ch_informationFlow() into per-cell median P-values. 
%
% Updated:
% J Noh, 2021/02/05. Rename and add comments. 
% Jungsik Noh, 2020/01/15


% load ML, example MD
ML.getMovies()
md1 = ML.getMovie(1);
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])
 
fsaveName0 = ['intraCh', num2str(iChan1)];
disp(fsaveName0)

%% input parsing

MDs = ML.getMovies();
num = numel(MDs);

ch1ActmapName = [chan1Name]; 
 
ip = inputParser;
ip.addParameter('outDirName', ['GC_SPAR1ch_', fsaveName0]);
%ip.addParameter('timeInterval', md1.timeInterval_);
%ip.addParameter('ccfLagMax', -1);
parse(ip, varargin{:})
p = ip.Results;

fname_GCresult = ['GC_', fsaveName0, '_winFits_Layers.mat'];
%MDtimeIntvl = p.timeInterval;

%% outDir

outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%%  GC winFits, winFits{direction}{indL}

MDs = ML.getMovies();
num = numel(MDs);
%
cellLabels = cell(num, 1);
FpArr = cell(num, 1);
medFpvalMat = nan(num, 4, maxLayer);
medFDRMat = nan(num, 4, maxLayer);
medpartR2Mat = nan(num, 4, maxLayer);

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
    wmaxi = size(winFits{1}{1}, 1);
    effWinInd = ones(wmaxi, 1);
    effWinInd(omittedWin) = 0;
    
    cellLabels{i} = cellName;
    %copyfile(fullfile(mdDir, analName, ['PW_xcCurve_', ch12fname, '.png']), ...
    %    fullfile(outDir, [cellLabels{i}, '_PW_xcorrCurve_', ch12fname, '.png']) )
    S2 = load(fullfile(mdDir, analName, [fsaveName0, '_Fpvec_FDR.mat']));
    FpArr{i} = S2.FpArr;
    
    for indL = 1:maxLayer
        for direction = 1:4
            Fpvec = S2.FpvecLyrDrctn{indL, direction};
            medFpvalMat(i, direction, indL) = median(Fpvec, 'omitnan');    % median P-values
            partR2vec = winFits{direction}{indL}(:, 11);
            medpartR2Mat(i, direction, indL) = median(partR2vec, 'omitnan'); % median partial R-square
            fdrbh = S2.FDRLyrDrctn{indL, direction};
            medFDRMat(i, direction, indL) = median(fdrbh, 'omitnan');
            
            effSize = sum(~isnan(Fpvec) & effWinInd);
            %disp(effSize)
            if effSize <  6                %0.5*sum(effWinInd)
                medFpvalMat(i, direction, indL) = nan;
                medpartR2Mat(i, direction, indL) = nan;
            end
        end
    end
    
end

log10medFpvalMat = log10(medFpvalMat);
%%
medFpvalMatfull = [];
for indL = 1:maxLayer
    medFpvalMatfull = [medFpvalMatfull, medFpvalMat(:, :, indL), nan(num, 1)];
end
medFDRMatfull = [];
for indL = 1:maxLayer
    medFDRMatfull = [medFDRMatfull, medFDRMat(:, :, indL), nan(num, 1)];
end
medpartR2Matfull = [];
for indL = 1:maxLayer
    medpartR2Matfull = [medpartR2Matfull, medpartR2Mat(:, :, indL), nan(num, 1)];
end


%
save(fullfile(outDir, ['GC_results_', fsaveName0, '.mat']), ...
    'cellLabels', 'FpArr', 'ch1ActmapName', ...
    'maxLayer', 'analName', 'medFpvalMat', 'medFDRMat', ...
    'log10medFpvalMat', 'medpartR2Mat', 'medFpvalMatfull', 'medFDRMatfull', 'medpartR2Matfull')

%% writetable medFpvalMat

directionName = {'Left', 'Right', 'Outside', 'Inside'};
a = 1:maxLayer;
%a1 = cellstr(num2str(a'));
a1 = string(a');
condNames = strcat(a1, 'L')';

directionName2 = {'Left', 'Right', 'Outside', 'Inside', ''};
condNames2 = repmat(directionName2, 1, maxLayer);

tab = array2table(medFpvalMatfull);
tab2 = array2table(medpartR2Matfull);
tab3 = array2table(medFDRMatfull);

writetable(tab, fullfile(outDir, ['GC_medFpval_Table_', fsaveName0, '.csv']))
writetable(tab2, fullfile(outDir, ['GC_medpartialR2_Table_', fsaveName0, '.csv']))
writetable(tab3, fullfile(outDir, ['GC_medFDR_Table_', fsaveName0, '.csv']))
%writetable(tab, 'upperQuartilePersTime.csv')

%%  signrank test p-val, 2019/05/16, nan checked

signRankPvec = nan(1, size(medFpvalMatfull,2));
for l = 1:size(medFpvalMatfull, 2)
    pvals = medFpvalMatfull(:, l);
    if any(~isnan(pvals));
        signRankPvec(l) = signrank(pvals, 0.05, 'tail', 'left');
    else
        signRankPvec(l) = NaN;
    end
end
signRankPvec1 = round(signRankPvec, 3);
medmedFpvalVec = round(median(medFpvalMatfull, 1, 'omitnan'), 3);

%% significance indicator

sigIndic = num2cell((signRankPvec1 < 0.05)');

%%  myBoxplot

fb = figure;
log10medFpvalMatfull = log10(medFpvalMatfull);
%matOut = max(-3, -log10medFpvalMat);
matOut = -log10medFpvalMatfull;
matOut = min(6, matOut, 'includenan');
%
if (size(matOut, 1) > 1)
    f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames2); %, 'LabelOrientation', 'inline');
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
ylabel('Median P-value (-log10)')

ax.XTickLabel = condNames2;
ax.XTickLabelRotation = 45;
ax.FontSize = 10;

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

title0 = ['GC per-cell median P-values, ', chan1Name];
%title1 = ['from ', chan1Name, ' to ', chan2Name, ' given ', chan3Name];
title2 = ['GC Pval: ', num2str(medmedFpvalVec)];
title3 = ['rankTest Pval: (', num2str(signRankPvec1), ')'];
title({title0;  title2; title3})

legend(s, cellLabels, 'Location', 'eastoutside')

%%
saveas(fb, fullfile(outDir, ['GC_BPFpval2_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['GC_BPFpval2_', fsaveName0, '.fig']), 'fig')


%% remove later....
%% rejection region of H1: mean(Median P-val) < 0.05 (2019/04, PvalLogScale) (removed)

%%  myBoxplot, medpartR2Mat

fb = figure('Visible', 'on');

% low bound 3
matOut = 100 * medpartR2Matfull;
ymaxR2 = max(10, max(matOut(:)));  % partR2 max is at least 10%

%
if (size(matOut, 1) > 1)
    f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames2);
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

ax.XTickLabelRotation = 45;
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

title0 = ['GC partial R-squares, ', chan1Name];
%
%title1 = num2str(medmedpartR2Mat);
title({title0; title1})

legend(s, cellLabels, 'Location', 'eastoutside')


%%
saveas(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['GC_BP_partialRsquare_', fsaveName0, '.fig']), 'fig')

%%  myBoxplot,  medFDRMatfull

fb = figure('Visible', 'on');

% low bound 3
matOut = medFDRMatfull;
%ymaxR2 = max(10, max(matOut(:)));  % partR2 max is at least 10%

%
if (size(matOut, 1) > 1)
    f1 = boxplot(matOut, 'Whisker', Inf, 'Labels', condNames2);
end
%ylim([0, max(matOut(:))*1.1])
ylim([0, 1])

%h = refline([0 -log10(0.05)]);h.Color='k'; h.LineStyle='--';
%h1 = refline([0 6]);  %h.Color='k'; h.LineStyle='--';
%ptick = [1 0.5 0.05, 0.01, 0.001 0.0001 0.000001];
%logptick = -log10(ptick);
ax = gca;
%ax.YTickMode = 'manual';
%ax.YTick = logptick;
%ax.YTickLabel = ptick;    %{'0.05', '0.01', '0.001', ' '};
ylabel('Median FDR')

ax.XTickLabelRotation = 45;
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
title1 = sprintf("%02.1f \t", medmedpartR2Mat(:));

title0 = ['GC FDR, ', chan1Name];
%
%title1 = num2str(medmedpartR2Mat);
title({title0; title1})

legend(s, cellLabels, 'Location', 'eastoutside')


%%
saveas(fb, fullfile(outDir, ['GC_BP_FDR_', fsaveName0, '.png']), 'png')
saveas(fb, fullfile(outDir, ['GC_BP_FDR_', fsaveName0, '.fig']), 'fig')

%%
disp('==== MLsummary_GC_SPAR1ch_informationFlow is finished!! ====')
disp('== :) close all ')

end
