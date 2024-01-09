function MLsummary_XcorrCurvesVelAcf(ML, iChan1, iChan2, chan1Name, chan2Name, ...
    maxLayer, analNameAcf, analNameXcf, varargin)
% MLsummary_XcorrCurvesVelAcf Collect/Summarize cross correlation curves
% for movies computed and saved by mapXcorrCurvePermutation.m. 
% It computed 95% confidence bands of the Xcorrelations at each lag using 2*SEM 
% under the assumption that movies are independent replicates.
%
% Usage:
%       MLsummary_XcorrCurvesVelAcf(ML, 1, 0, 'Actin', 'Vel', 3, ...
%           'mapDescriptives', 'mapCrossCorr')
%
% Input:
%       ML          - a movieList object
%       iChan1      - the 1st channel index
%       chan1Name   - a short name for channel1.
%       iChan2      - the 2nd channel index
%       chan2Name   - a short name for channel2.
%       maxLayer    - maximum layer to be analyzed 
%       analNameAcf - the folder name for output for edge velocity
%                       (iChan=0) to collect acf curves
%       analNameXcf - the folder name for output from
%                     mapXcorrCurvePermutation.m to collect xcf curves
%
% Output: .png/.fig/.mat files are saved in the ML.outputDirectory_/acfXcf_Ch1Ch0 by default. 
%       
% Option:
%       outDirName  - Specify a name of output directory.
%
% Updates:
% J Noh, 2021/07/05. Adjust lagSizeAcf (Xcf) in case of 'data driven'
% option (p.lagMax0 == 5).
% J Noh, 2018/11/14. Edit legend for layers.
% J Noh, 2018/10/24. Change lagMax0 default from nan to 0 to be consistent
% with GUI. Add a 'figFlag' option. 
% J Noh, 2018/04/03. Debug all NaN case.
% J Noh, 2018/01/23. Main plots are now saved as three formats.
% J Noh, 2017/11/05. 
% The ouput name of autocorr for Vel is now 'Chan0_Avg_autocorLayers.mat', 
% which is in general format from mapDescriptives_OneChan(iChan=0).
% J Noh, 2017/09/25. Include the summary of ACF of channels.
% Jungsik Noh, 2017/05/17
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
%
% This file is part of GrangerCausalityAnalysisPackage.
% 
% GrangerCausalityAnalysisPackage is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GrangerCausalityAnalysisPackage is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GrangerCausalityAnalysisPackage.  If not, see <http://www.gnu.org/licenses/>.
% 
% 


disp('===============================================================')
disp(['== Cross correlations between ', chan1Name, ' and ', chan2Name])
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
ip.addParameter('outDirName', ['acfXcf_', ch1fname, ch2fname]);
ip.addParameter('timeInterval', md1.timeInterval_);
ip.addParameter('lagMax0', 5);
ip.addParameter('figFlag', 'off');


parse(ip, varargin{:})
p = ip.Results;


fname_xcorrMat = ['xcorrMat_', ch1fname, ch2fname, '.mat'];
MDtimeIntvl = p.timeInterval;

%
set(groot,'defaultLegendAutoUpdate','off') 
 
%% setting up parameters (1)

outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%% get Avg_autocorLayers, cellNames, and copy cell ACF plots

acfvecsize = nan(num, 1);
cellLabels = cell(num, 1);
acfCell = cell(num, 1);

for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    
    [folderName, cellLab0, ~] = fileparts(mdDir);
    cellName = [cellLab0(1:end)];    % Generic cellName
    %%%%    
    cellLabels{i} = cellName;
    copyfile(fullfile(mdDir, analNameAcf, 'acCurve_Chan0_1L.png'), ...
        fullfile(outDir, [cellLabels{i}, '_acCurve_Chan0_1L.png']) )    

    %%%% input
    load(fullfile(mdDir, analNameAcf, 'Chan0_Avg_autocorLayers.mat'));  % Avg_autocor    
    
    acfCell{i} = Avg_autocorLayers{1};
    acfvecsize(i) = size(Avg_autocorLayers{1}, 2);
end

acfvecsize0 = min(acfvecsize);
acfArr2 = nan(num, acfvecsize0);
for i=1:num
    acfArr2(i,:) = acfCell{i}(1:acfvecsize0);
end

MLmeanACF = mean(acfArr2, 1, 'omitnan');
[~, minid] = min(MLmeanACF);
tLag = [0:acfvecsize0-1].* MDtimeIntvl;
acfminTime = tLag(minid);
estProtRetCycleFrs = 2*(minid-1);

%% determin lagSizeAcf, lagSizeXcf

if ~(p.lagMax0 == 5)
    lagSizeAcf = 2*p.lagMax0 + 1;
else
    % lagSizeAcf = min(acfvecsize);     % old
    lagSizeAcf = min( 2 * estProtRetCycleFrs + 1, min(acfvecsize));
end

disp(['Lag size (frames) of ACF: ', num2str(lagSizeAcf-1), ' frames'])

%
xcfmatsize = nan(num, 1);
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    %%%% input
    load(fullfile(mdDir, analNameXcf, fname_xcorrMat));  % xcorrMat
    xcfmatsize(i) = size(xcorrMat{1}, 2);
end

% if ~isnan(p.lagMax0)
if ~(p.lagMax0 == 5)
    lagSizeXcf = 2*p.lagMax0 + 1;
else
    lagSizeXcf = min( 2 * estProtRetCycleFrs + 1, min(xcfmatsize));
end
 
disp(['Length (frames) of XCorrCurves: ', num2str(lagSizeXcf), ' frames'])


%%  Acf_vel

acfArr = nan(num, lagSizeAcf);

for i = 1:num
   
        xcmean = acfCell{i};
        %tmplen = numel(xcmean);
        if (numel(xcmean) < lagSizeAcf)
            xcmeanMiddle = [xcmean, nan(1, lagSizeAcf - numel(xcmean))];
        else
            xcmeanMiddle = xcmean(1:lagSizeAcf);
        end
            
        acfArr(i, :) = reshape(xcmeanMiddle, 1, []);
    
end

totavgXcorrCurve = squeeze(mean(acfArr, 1, 'omitnan')); 
 

%%
tLag = [0:lagSizeAcf-1].* MDtimeIntvl;

    fAcf = figure('Visible', p.figFlag);   %%  name
    
    acmeanPerCellindL = squeeze(acfArr(:, :))';
    
    p1 = plot(tLag, totavgXcorrCurve, 'r');
    p1.LineWidth = 2;
    hold on
    
    plot(tLag, acmeanPerCellindL, 'Color', 'k', 'LineWidth', 0.5)
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    
    %title('Vel')
    title(['Vel (ACF): meanACF min at ', num2str(acfminTime), ' s'] )
    
    xlabel('Time lag (s)');ylabel('Correlation')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    ax.FontSize = 14;
    
%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
ylim([-0.3, 0.3])
xlim([min(tLag), max(tLag)])


%% saveas

saveas(fAcf, fullfile(outDir, ['/acfVel', '.png']), 'png')
saveas(fAcf, fullfile(outDir, ['/acfVel', '.fig']), 'fig')


%%  Acf_ChanX

for iChan = [iChan1, iChan2]
    if (iChan == iChan1)
        chanName = chan1Name;
    else
        chanName = chan2Name;
    end
    
    if (iChan ~= 0)
        fname0 = ['Chan', num2str(iChan)];
        
MDs = ML.getMovies();
num = numel(MDs);
%num = num    %%%  input
cellLabels = cell(num, 1);

acfArr = cell(maxLayer, 1);
for indL = 1:maxLayer
    acfArr{indL} = nan(num, lagSizeAcf);
end


%
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    %cellName = [cellLab0(1:14)];   %% Adjust name extraction
    cellName = [cellLab0(1:end)];    % Generic cellName
    %[~, cellName] = fileparts(folderName);

    %%%% input
    load(fullfile(mdDir, analNameAcf, [fname0, '_Avg_autocorLayers.mat']));   

    %%%%
    
    %cellLabels{i} = cellName;
    %copyfile(fullfile(mdDir, analNameAcf, 'acCurveChan0.png'), ...
    %    fullfile(outDir, [cellLabels{i}, '_acCurveChan0.png']) )    
    
    %xcorrMat = Avg_autocor;
        %xcmean = mean(xcorrMat_tmp{indL}, 1, 'omitnan');
    
    xcmean = [];
    for indL = 1:maxLayer        
        xcmean = Avg_autocorLayers{indL};
        %tmplen = numel(xcmean);
    
        if (numel(xcmean) < lagSizeAcf)
            xcmeanMiddle = [xcmean, nan(1, lagSizeAcf - numel(xcmean))];
        else
            xcmeanMiddle = xcmean(1:lagSizeAcf);
        end
            
        acfArr{indL}(i, :) = reshape(xcmeanMiddle, 1, []);
    end
        
end

totavgXcorrCurve = cell(maxLayer, 1);
for indL = 1:maxLayer
totavgXcorrCurve{indL} = squeeze(mean(acfArr{indL}, 1, 'omitnan'));
end

 
 

%%
tLag = [0:lagSizeAcf-1].* MDtimeIntvl;

fAcf = cell(maxLayer, 1);
for indL = 1:maxLayer
    
    fAcf{indL} = figure('Visible', p.figFlag);   %%  name
    
    acmeanPerCellindL = squeeze(acfArr{indL}(:, :))';
    
    p1 = plot(tLag, totavgXcorrCurve{indL}, 'r');
    p1.LineWidth = 2;    
    hold on
    
    plot(tLag, acmeanPerCellindL, 'Color', 'k', 'LineWidth', 0.5)
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    
    title([chanName, '-', num2str(indL), 'L (ACF)'] )
    
    xlabel('Time lag (s)');ylabel('Correlation')
    set(gca, 'XGrid', 'on')
    ax = gca;
    ax.FontSize = 14;
    
%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
ylim([-1, 1])
xlim([min(tLag), max(tLag)])



%% saveas

saveas(fAcf{indL}, fullfile(outDir, [fname0, '_acf_', num2str(indL), 'L.png']), 'png')
saveas(fAcf{indL}, fullfile(outDir, [fname0, '_acf_', num2str(indL), 'L.fig']), 'fig')

end

    end
end


%%  Xcf

MDs = ML.getMovies();
num = numel(MDs);

%
num = num   %%%  input

cellLabels = cell(num, 1);

xcorrArr = nan(num, lagSizeXcf, maxLayer);
 
 
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
    copyfile(fullfile(mdDir, analNameXcf, ['xcCurve_', ch12fname, '_permT.png']), ...
        fullfile(outDir, [cellLabels{i}, '_xcorrCurve_', ch12fname, '_permT.png']) )    
    
    for indL = 1:maxLayer
        %xcmean = mean(xcorrMat{indL}, 1, 'omitnan');
        if all(isnan(xcorrMat{indL}))
            xcmean = nan(1, size(xcorrMat{indL}, 2));
        else
            xcmean = smoothingSplineCorMap(xcorrMat{indL});
        end
            
        tmplen = numel(xcmean);
        if (tmplen < lagSizeXcf)
            disp(['== Lag size of xcorrMat is too short in movie #', num2str(i), '. =='])
        end
        xcmeanMiddle = xcmean(1+(tmplen - lagSizeXcf)/2:lagSizeXcf+(tmplen - lagSizeXcf)/2);
        xcorrArr(i, :, indL) = reshape(xcmeanMiddle, [], 1);
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

%colOrd = get(gca,'ColorOrder');
%colOrd1 = colOrd(1:7, :);                   % Revise later together with legend.

f1 = cell(maxLayer,1);
tmpYlim = nan(maxLayer, 2);
% 2018/11. For GUI.
min0 = min(min(xcorrArr(:)), 0) * 1.1;
max0 = max(max(xcorrArr(:)), 0) * 1.1;

for indL = 1:maxLayer

    f1{indL} = figure('Visible', p.figFlag);   %%  name
    
    xcmeanPerCellindL = squeeze(xcorrArr(:, :, indL))';
    
    p1 = plot(1:lagSizeXcf, totavgXcorrCurve(:, indL), 'r');
    p1.LineWidth = 2;
    
    hold on
    
    plot(1:lagSizeXcf, xcmeanPerCellindL, 'Color', 'k', 'LineWidth', 0.5)
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, '_t) -', num2str(indL), 'L'])
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
%    figure(f1{indL});
%    ylim([min0, max0])
%end




%% saveas
for indL=1:maxLayer
saveas(f1{indL}, fullfile(outDir, ['/xcorrSummary_', num2str(indL), 'L.png']), 'png')
saveas(f1{indL}, fullfile(outDir, ['/xcorrSummary_', num2str(indL), 'L.fig']), 'fig')
end

for indL=1:maxLayer
    figure(f1{indL});
    title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)
    
    saveas3format(f1{indL}, outDir, [ch1ActmapName, '_xcorrSummary_', num2str(indL), 'L']);
    

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

f2 = cell(maxLayer,1);
tmpYlim = nan(maxLayer, 2);
% 2018/11. For GUI.
min0 = min(min(xcorrArr(:)), 0) * 1.1;
max0 = max(max(xcorrArr(:)), 0) * 1.1;

for indL = 1:maxLayer

    f2{indL} = figure('Visible', p.figFlag);   %%  name
    
    xcMVecsIndL = squeeze(xcorrArr(:, :, indL));
    xcTotM = mean(xcMVecsIndL, 1, 'omitnan');
    sem = std(xcMVecsIndL, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(1:lagSizeXcf, xcTotM, 2*sem, 'lineprops', '-r');
    s1.mainLine.LineWidth = 2;
    
    hold on
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, '_t) -', num2str(indL), 'L'])
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
saveas(f2{indL}, fullfile(outDir, ['/xcorrMean_', num2str(indL), 'L.png']), 'png')
saveas(f2{indL}, fullfile(outDir, ['/xcorrMean_', num2str(indL), 'L.fig']), 'fig')
end

for indL=1:maxLayer
    figure(f2{indL});
    title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)
    
    saveas3format(f2{indL}, outDir, [ch1ActmapName, '_xcorrMean_', num2str(indL), 'L']);

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
CorMaxLags = nan(num, 2*maxLayer);
%rownames = {'cm1L', 'lag1L', 'cm2L', 'lag2L', 'cm3L', 'lag3L'};
%rownames = {'CM1L', 'lag1L'};
rownames = cell(1, 2*maxLayer);
for indL = 1:maxLayer
    rownames{1 + 2*(indL-1)} = ['CM', num2str(indL), 'L'];
    rownames{2*indL} = ['lag', num2str(indL), 'L'];
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
    
    cmCell = xcorrMatToCMLag(xcorrMat);
    cmMat = cell2mat(cmCell)';
    cmVec = reshape(cmMat, 1, []);
    
    for j = 1:numel(cmVec)
        CorMaxLags(i, j) = cmVec(j);
    end
    cellLabels{i} = cellName;
end

T0 = table(cellLabels);
T1 = array2table(CorMaxLags(:, 1:2*maxLayer));
T1.Properties.VariableNames = rownames;
T = [T0, T1];

%%
save(fullfile(outDir, 'CorrMaximaLag_table.mat'), 'T', 'cellLabels', 'CorMaxLags')



%% cor maxima/minima plot

fcm = cell(maxLayer, 1);

ind0 = 1:maxLayer;
cors = CorMaxLags(:, 1+2*(ind0-1));
lags = CorMaxLags(:, 2*ind0);

corsMin = min(min(cors(:)), -0.01) * 1.2;
corsMax = max(max(cors(:)), 0.01) * 1.2;
lagsMin = min(lags(:)) - 3;
lagsMax = max(lags(:)) + 3;


for indL=1:maxLayer

    r1 = CorMaxLags(:, 1+2*(indL-1));

    l1 = CorMaxLags(:, 2*indL) .*MDtimeIntvl;

    % xlim0 = [-15, 50]
    % ylim([-0.1 0.15]) 

    fcm{indL} = figure('Visible', p.figFlag);   %%  name
    scatter(l1, r1, [], 'r', 'filled')
    text(l1+0.5, r1+rand(num,1)*0.02-0.01, cellLabels)
    ylabel('Corr Maximum/Minmum')
    xlabel('Lag (relative to vel) (sec)')
    
    rMean = mean(r1);
    lMean = mean(l1);
    lMedian = median(l1);
    rMedian = median(r1);
    title1 = ['Layer ', num2str(indL)];
    title2 = ['avgLag: ', num2str(lMean), '(', num2str(lMedian), ')', ' avgCor: ', num2str(rMean), '(', num2str(rMedian), ')'];
    title({title1, title2})
    hold on
    scatter(lMean, rMean, [], '+')

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
saveas(fcm{indL}, fullfile(outDir, ['/CMLag', num2str(indL), 'L.png']), 'png')
saveas(fcm{indL}, fullfile(outDir, ['/CMLag', num2str(indL), 'L.fig']), 'fig')
end



%%  xcorr summary2

lagMax = (lagSizeXcf - 1)/2;
lagGrid = floor(lagMax/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeIntvl, 2);

%colOrd = get(gca,'ColorOrder');
%colOrd1 = colOrd(1:7, :);                   % Revise later together with legend.

f1 = cell(maxLayer,1);
tmpYlim = nan(maxLayer, 2);
% 2018/11. For GUI.
min0 = min(min(xcorrArr(:)), 0) * 1.1;
max0 = max(max(xcorrArr(:)), 0) * 1.1;

for indL = 1:maxLayer

    f1{indL} = figure('Visible', p.figFlag);   %%  name
    
    xcmeanPerCellindL = squeeze(xcorrArr(:, :, indL))';
    

    hold on
    p1 = plot(1:lagSizeXcf, totavgXcorrCurve(:, indL), 'r');
    p1.LineWidth = 2;
    
    xcMVecsIndL = squeeze(xcorrArr(:, :, indL));
    xcTotM = mean(xcMVecsIndL, 1, 'omitnan');
    sem = std(xcMVecsIndL, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(1:lagSizeXcf, xcTotM, 2*sem, 'lineprops', '-r');
    s1.mainLine.LineWidth = 2;
    s1.edge(1).LineWidth = 0.5; s1.edge(2).LineWidth = 0.5;
    
    pp = plot(1:lagSizeXcf, xcmeanPerCellindL, 'Color', 'k', 'LineWidth', 0.5);
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, '_t) -', num2str(indL), 'L'])
    xlabel('Time lag (s)');ylabel('Correlation')

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
%    figure(f1{indL});
%    ylim([min0, max0])
%end

 
%% saveas 2
for indL=1:maxLayer
saveas(f1{indL}, fullfile(outDir, ['xcorrSummary2_', num2str(indL), 'L.png']), 'png')
saveas(f1{indL}, fullfile(outDir, ['xcorrSummary2_', num2str(indL), 'L.fig']), 'fig')
end

for indL=1:maxLayer
    figure(f1{indL});
    title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)
    
    saveas3format(f1{indL}, outDir, [ch1ActmapName, '_xcorrSummary2_', num2str(indL), 'L'])

end
if strcmp(p.figFlag, 'off')
for indL=1:maxLayer
    close(f1{indL})
end
end

 

%%
disp('==== MLsummary_XcorrCurvesVelAcf is finished!! ====')

end
