function MLsummary_layerCrossCorrCurves(ML, iChan1, iChan2, chan1Name, chan2Name, ...
    maxLayer, analNameXcf, varargin)
% MLsummary_layerCrossCorrCurves Collect/Summarize inter-layer cross correlation curves 
% for movies computed and saved by layerCrossCorrCurvePermutation.m. 
% It computed 95% confidence bands of the Xcorrelations at each lag using 2*SEM 
% under the assumption that movies are independent replicates.
%
% Usage:
%       MLsummary_layerCrossCorrCurves(ML, 1, 0, 'Actin', 'Vel', 3, ...
%           'mapDescriptives', 'mapCrossCorr')
%
% Input:
%       ML          - a movieList object
%       iChan1      - the 1st channel index
%       chan1Name   - a short name for channel1.
%       iChan2      - the 2nd channel index
%       chan2Name   - a short name for channel2.
%       maxLayer    - maximum layer to be analyzed 
%       analNameXcf - the folder name for output from
%                     fieldXcorrCurvePermutation.m to collect xcf curves
%
% Output: .png/.fig/.mat files are saved in the ML.outputDirectory_/FieldXcf_Ch1Ch0 by default. 
%       
% Option:
%       outDirName  - Specify a name of output directory.
%
% Updates:
% J Noh, 2021/02/05. Rename and add comments. 
% Jungsik Noh, 2019/03/26


disp('===============================================================')
disp(['== Inter-layer cross correlations between ', chan1Name, ' and ', chan2Name])
disp('===============================================================')

% load ML, example MD
ML.getMovies();
md1 = ML.getMovie(1);
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])

%%
MDs = ML.getMovies();
num = numel(MDs);

ch1ActmapName = chan1Name;
ch2ActmapName = chan2Name;

ch1fname = ['Ch', num2str(iChan1)];
ch2fname = ['Ch', num2str(iChan2)];
ch12fname = [ch1fname, ch2fname];

ip = inputParser; 
ip.addParameter('outDirName', ['FieldXcf_', ch1fname, ch2fname]);
ip.addParameter('timeInterval', md1.timeInterval_);
ip.addParameter('lagMax0', 5);
ip.addParameter('figFlag', 'off');


parse(ip, varargin{:})
p = ip.Results;


fname_xcorrMat = ['xcorrMatArr_', ch1fname, ch2fname, '.mat'];
MDtimeIntvl = p.timeInterval;

%
set(groot,'defaultLegendAutoUpdate','off') 
 
%% setting up parameters

outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%acfvecsize = nan(num, 1); ...

%
xcfmatsize = nan(num, 1);
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    %%%% input
    load(fullfile(mdDir, analNameXcf, fname_xcorrMat));  % xcorrMatArr
    xcfmatsize(i) = size(xcorrMatArr{1,1}, 2);
end

if ~isnan(p.lagMax0)
    lagSizeXcf = 2*p.lagMax0 + 1;
else
    lagSizeXcf = min(xcfmatsize);
end

% or manually assign a value like: lagSizeXcf = 81;
disp(['Length (frames) of XCorrCurves: ', num2str(lagSizeXcf), ' frames'])



%%  Acf_vel


%%  Acf_ChanX


%%  Xcf

MDs = ML.getMovies();
num = numel(MDs);

%
num = num   %%%  input

cellLabels = cell(num, 1);

xcorrArr = nan(num, lagSizeXcf, maxLayer, maxLayer);
 
 
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    cellName = [cellLab0(1:end)];    % Generic cellName
    %[~, cellName] = fileparts(folderName);

    %%%% input
    load(fullfile(mdDir, analNameXcf, fname_xcorrMat));   % xcorrMat_tmp
    %xcorrMat = xcorrMat_tmp;
    %%%%
    %pooledPermXcorrMean{i} = nan([size(permXcorrMean{1}), maxLayer]); 
    
    cellLabels{i} = cellName;
    
    for indL2 = 1:maxLayer
    copyfile(fullfile(mdDir, analNameXcf, ['xcCurve_', ch12fname, '_to', num2str(indL2), 'L.png']), ...
        fullfile(outDir, [cellLabels{i}, '_xcorrCurve_', ch12fname, '_to', num2str(indL2), 'L.png']) ) 
    for indL = 1:maxLayer    
        %xcmean = mean(xcorrMat{indL}, 1, 'omitnan');
        if all(isnan(xcorrMatArr{indL, indL2}))
            xcmean = nan(1, size(xcorrMatArr{indL, indL2}, 2));
        else
            %xcmean = smoothingSplineCorMap(xcorrMatArr{indL, indL2});
            xcmean = nanmean(xcorrMatArr{indL, indL2}, 1);
        end
            
        tmplen = numel(xcmean);
        if (tmplen < lagSizeXcf)
            disp(['== Lag size of xcorrMat is too short in movie #', num2str(i), '. =='])
        end
        xcmeanMiddle = xcmean(1+(tmplen - lagSizeXcf)/2:lagSizeXcf+(tmplen - lagSizeXcf)/2);
        xcorrArr(i, :, indL, indL2) = reshape(xcmeanMiddle, [], 1);
    end
    end

end


totavgXcorrCurve = squeeze(mean(xcorrArr, 1, 'omitnan'));


%
save(fullfile(outDir, 'aggregatedXcorr.mat'), ...
      'cellLabels', 'xcorrArr', 'ch1ActmapName', ...
      'ch2ActmapName', 'maxLayer', 'analNameXcf', 'lagSizeXcf')



%% 

lagMax = (lagSizeXcf - 1)/2;
lagGrid = floor(lagMax/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeIntvl, 2);

ftmp = figure;
colOrd = get(gca,'ColorOrder');
close(ftmp)

nRow0 = size(colOrd,1); multipleNum = ceil(maxLayer/nRow0);
colOrdExt = repmat(colOrd, multipleNum, 1);
colOrd1 = colOrdExt(1:maxLayer, :);  

f1 = cell(maxLayer,1);
tmpYlim = nan(maxLayer, 2);
plobj = cell(maxLayer, 1);
% 2018/11. For GUI.
min0 = min(min(xcorrArr(:)), 0) * 1.1;
max0 = max(max(xcorrArr(:)), 0) * 1.1;

for indL2 = 1:maxLayer

    f1{indL2} = figure('Visible', p.figFlag);   %%  name

    hold on
    
for indL = 1:maxLayer

    figure(f1{indL2});
    xcmeanPerCellindL = squeeze(xcorrArr(:, :, indL, indL2))';
    plot(1:lagSizeXcf, xcmeanPerCellindL, 'Color', colOrd1(indL,:), 'LineWidth', 0.5)
    
    plobj{indL} = plot(1:lagSizeXcf, totavgXcorrCurve(:, indL, indL2), 'Color', colOrd1(indL,:), ...
        'DisplayName', [num2str(indL), 'L']);
    plobj{indL}.LineWidth = 2;
    
    %hold on
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
end

    figure(f1{indL2});
    legend([plobj{:}], 'Location','northoutside','Orientation','horizontal')

    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, num2str(indL2),  'L_t)'])
    
    xlabel('Time lag (s)');ylabel('Correlation')
%    legend('1L', '2L', '3L', '4L',  'Location','northoutside','Orientation','horizontal');
    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    %h = refline([0 0]); h.Color = [.5 .5 .5];
    %h = refline([Inf 0]); h.Color = [.5 .5 .5];
    set(gca, 'XGrid', 'on')
    
ax = gca;
ax.FontSize = 14;

xlim([1, lagSizeXcf])
    %tmpYlim(indL, :) = ax.YLim;
    ylim([min0, max0])
    pause(1)  
end

%min0 = min(tmpYlim(:, 1));
%max0 = max(tmpYlim(:, 2));
%for indL = 1:maxLayer
%    figure(f1{indL});
%    ylim([min0, max0])
%end




%% saveas
for indL=1:maxLayer
saveas(f1{indL}, fullfile(outDir, ['fieldXcorrSummary_to', num2str(indL), 'L.png']), 'png')
saveas(f1{indL}, fullfile(outDir, ['fieldXcorrSummary_to', num2str(indL), 'L.fig']), 'fig')
end

for indL=1:maxLayer
    figure(f1{indL});
    title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)
    
    saveas3format(f1{indL}, outDir, [ch1ActmapName, '_fieldXcorrSummary_to', num2str(indL), 'L']);
    

    pause(0.1)
end
if strcmp(p.figFlag, 'off')
for indL=1:maxLayer
    close(f1{indL})
end
end


 
%%  integrate movies SEM

lagMax = (lagSizeXcf - 1)/2;
lagGrid = floor(lagMax/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeIntvl, 2);

%colOrd = get(gca,'ColorOrder');
%colOrd1 = colOrd(1:5, :);                   % Revise later together with legend.
ftmp = figure;
colOrd = get(gca,'ColorOrder');
close(ftmp)

nRow0 = size(colOrd,1); multipleNum = ceil(maxLayer/nRow0);
colOrdExt = repmat(colOrd, multipleNum, 1);
colOrd1 = colOrdExt(1:maxLayer, :);  

f2 = cell(maxLayer,1);
tmpYlim = nan(maxLayer, 2);
plobj = cell(maxLayer, 1);
% 2018/11. For GUI.
min0 = min(min(xcorrArr(:)), 0) * 1.1;
max0 = max(max(xcorrArr(:)), 0) * 1.1;

for indL2 = 1:maxLayer

    f2{indL2} = figure('Visible', p.figFlag);   %%  name

for indL = 1:maxLayer
    
    xcMVecsIndL = squeeze(xcorrArr(:, :, indL, indL2));
    xcTotM = mean(xcMVecsIndL, 1, 'omitnan');
    sem = std(xcMVecsIndL, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(1:lagSizeXcf, xcTotM, 2*sem, 'lineprops', {'Color', colOrd1(indL,:)}); %, 'markerfacecolor', colOrd1(indL2,:)});
    plobj{indL} = s1.mainLine;
    plobj{indL}.LineWidth = 2;
    plobj{indL}.DisplayName = [num2str(indL), 'L'];
    %s1.mainLine.Color = colOrd1(indL2,:);
    hold on
end

    %a = 1:maxLayer; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'L'});
    %legend(d, 'Location','northoutside','Orientation','horizontal');
    legend([plobj{:}],  'Location','northoutside','Orientation','horizontal')
    
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, num2str(indL2),  'L_t)'])
    
    xlabel('Time lag (s)');ylabel('Correlation')
%    legend('1L', '2L', '3L', '4L',  'Location','northoutside','Orientation','horizontal');
    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    %h = refline([0 0]); h.Color = [.5 .5 .5];
    %h = refline([Inf 0]); h.Color = [.5 .5 .5];
    set(gca, 'XGrid', 'on')
    
ax = gca;
ax.FontSize = 14;

xlim([1, lagSizeXcf])

    %tmpYlim(indL, :) = ax.YLim;
    ylim([min0, max0])

end

%min0 = min(tmpYlim(:, 1));
%max0 = max(tmpYlim(:, 2));
%for indL = 1:maxLayer
%    figure(f2{indL});
%    ylim([min0, max0])
%end



%% saveas
for indL=1:maxLayer
saveas(f2{indL}, fullfile(outDir, ['fieldXcorrMean_to', num2str(indL), 'L.png']), 'png')
saveas(f2{indL}, fullfile(outDir, ['fieldXcorrMean_to', num2str(indL), 'L.fig']), 'fig')
end

for indL=1:maxLayer
    figure(f2{indL});
    title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)
    
    saveas3format(f2{indL}, outDir, [ch1ActmapName, '_fieldXcorrMean_to', num2str(indL), 'L']);

end
if strcmp(p.figFlag, 'off')
for indL=1:maxLayer
    close(f2{indL})
end
end




%%  CMLags

MDs = ML.getMovies();
num = numel(MDs);


%%

cellLabels = cell(num, 1);
%CorMaxLags = nan(num, 2*maxLayer);
CorMaxima = nan(num, maxLayer*maxLayer);
CorMaximaLags = nan(num, maxLayer*maxLayer);
%rownames = {'cm1L', 'lag1L', 'cm2L', 'lag2L', 'cm3L', 'lag3L'};
%rownames = {'CM1L', 'lag1L'};
rownames1 = cell(1, maxLayer*maxLayer);
rownames2 = cell(1, maxLayer*maxLayer);
for indL2 = 1:maxLayer
    for indL = 1:maxLayer
        rownames1{indL + maxLayer*(indL2-1)} = ['CM', num2str(indL), 'Lto', num2str(indL2)];
        rownames2{indL + maxLayer*(indL2-1)} = ['lags', num2str(indL), 'Lto', num2str(indL2)];
    end
    
end


for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    cellName = [cellLab0(1:end)];    % Generic cellName
    %[~, cellName] = fileparts(folderName);

    %%%% input
    %analName = 'mapCrossCorr_0207_Imp1L3';
    load(fullfile(mdDir, analNameXcf, fname_xcorrMat));
    %%%%
    
    cmCell = xcorrMatToCMLag_nmean(xcorrMatArr);
    cmMat = cell2mat(cmCell)';
    %cmVec = reshape(cmMat, 1, []);
    
    for j = 1:size(cmMat,2)
        CorMaxima(i, j) = cmMat(1,j);
        CorMaximaLags(i, j) = cmMat(2,j);
    end
    cellLabels{i} = cellName;
end

T0 = table(cellLabels);
%T1 = array2table(CorMaxLags(:, 1:2*maxLayer));
T1 = array2table(CorMaxima(:, :));
T1.Properties.VariableNames = rownames1;
T2 = array2table(CorMaximaLags(:, :));
T2.Properties.VariableNames = rownames2;
T = [T0, T1, T2];

%%
save(fullfile(outDir, 'CorrMaximaLag_table.mat'), 'T', 'cellLabels', 'CorMaxima', 'CorMaximaLags')



%% cor maxima/minima plot

fcm = cell(maxLayer, 1);
plobj = cell(maxLayer, 1);

ind0 = 1:maxLayer;
cors = CorMaxima(:, :);
lags = CorMaximaLags(:, :);

corsMin = min(min(cors(:)), -0.01) * 1.2;
corsMax = max(max(cors(:)), 0.01) * 1.2;
lagsMin = min(lags(:)) - 3;
lagsMax = max(lags(:)) + 3;


for indL2=1:maxLayer

    r1 = CorMaxima(:, ind0+maxLayer*(indL2-1));

    l1 = CorMaximaLags(:, ind0+maxLayer*(indL2-1)) .*MDtimeIntvl;

    % xlim0 = [-15, 50]
    % ylim([-0.1 0.15]) 

    fcm{indL2} = figure('Visible', p.figFlag);   %%  name
    
for indL = 1:maxLayer
    l10 = l1(:,indL); r10 = r1(:, indL);
    scatter(l10, r10, [], colOrd1(indL,:), 'filled')
    text(l10+0.5, r10+rand(num,1)*0.02-0.01, cellLabels)
    hold on
    
    rMean = mean(r10);
    lMean = mean(l10);
    lMedian = median(l10);
    rMedian = median(r10);
    plobj{indL} = scatter(lMean, rMean, [], '+');
    plobj{indL}.DisplayName = [num2str(indL), 'L'];
end

    ylabel('Corr Maximum/Minmum')
    xlabel('Lag (sec)')    
    legend([plobj{:}],  'Location','northoutside','Orientation','horizontal')

    title(['toLayer ', num2str(indL2)]);
    %title2 = ['avgLag: ', num2str(lMean), '(', num2str(lMedian), ')', ' avgCor: ', num2str(rMean), '(', num2str(rMedian), ')'];
    %title({title1, title2})

    ylim([corsMin, corsMax])
    xlim0 = [lagsMin lagsMax].* MDtimeIntvl;
    xlim(xlim0)

    h=refline([0,0]); %h.Color = [.5 .5 .5];
    h = line([0, 0], ylim); %h.Color = [.5 .5 .5];
    
    ax = gca;
    ax.FontSize = 14;

end

%% saveas
for indL=1:maxLayer
saveas(fcm{indL}, fullfile(outDir, ['/CMLag_to', num2str(indL), 'L.png']), 'png')
saveas(fcm{indL}, fullfile(outDir, ['/CMLag_to', num2str(indL), 'L.fig']), 'fig')
end

%%
disp('==== MLsummary_layerCrossCorrCurves is finished!! ====')
disp('== :) close all ')

end
