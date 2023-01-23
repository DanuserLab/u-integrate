function mapDescriptives_Vel_LB(MD, figuresDir, varargin)
% mapDescriptives_Vel_LB PERFORM Ljung-Box test (lbqtest.m) to edge velocity time
% series for each window. The test identifies whether the velocity TS is
% white-noise or not so that windows are classified into quiescent ones
% (when the velocity is a white-noise) and active ones (when the
% autocorrelations of the velocity are significant).
%
% Input:
%       MD          - movieData object
%       figuresDir  - a directory where plots are saved as png files
%
% Output:           - 'indActive_windowIndex.mat'. 
%                   - png files are saved in the figuresDir
%
% Option:
%       figFlag     - if 'on', matlab figures are ploted. Default is 'off'.
%       impute      - if true, moderate missing values are imputed by using
%                   knnimpute.m function.
%       omittedWindows  
%                   - window index in which activities will be replaced by
%                   NaN. Default is null.
%       subFrames
%                   - specified frames will be only used.  
%       derivative  - if true, it computes differenced velocities per
%                   second.
%       topograph   - if 'off' topographs are not plotted. Default is 'on'.
%
% Updated:
% J Noh, 2019/06/11. 'EWMA' is added. mapOutlierImputation() options
% updated.
% J Noh, 2019/03/25. checkWindowJump added.
% J Noh, 2019/02/27. Output 'Vel_meanACF_corAvg_active.mat' for corAvg_active. 
% J Noh, 2018/11/14. Edit legend for layers.
% J Noh, 2018/02/21. 'movmeanNum' option can control the number of
% neighboring P-values of LB tests. Default is 3. 
%
% Jungsik Noh, 2018/01/29.
%
% Copyright (C) 2023, Danuser Lab - UTSouthwestern 
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
 


ip = inputParser; 
ip.addParameter('figFlag', 'off');
ip.addParameter('impute', true);    
ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('derivative', false);
ip.addParameter('topograph', 'on');
ip.addParameter('movingAvgSmoothing', true);  
ip.addParameter('smParam', 0.9);
ip.addParameter('movmeanNum', 3);           % odd movmeanNum is preferred.
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing. 
ip.addParameter('outlSigma', 5);

ip.parse(varargin{:});
p = ip.Results;

%figFlag = p.figFlag;
set(groot,'defaultLegendAutoUpdate','off') 


%%  figuresDir setup
if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['mapDescriptives_Vel_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')

%%  getting Maps from channels

switch p.derivative
    case false
        chan0Title = 'Velocity (nm/sec)';
        chan0Name = 'Vel';
        iChan = 0;
    case true
        chan0Title = 'D.Velocity (nm/sec^2)';
        chan0Name = 'D.Vel';
        iChan = 10;
end

%  fnameChan = double(p.derivative)*10 + iChan;

        disp(chan0Name)
        disp(chan0Title)
        disp(iChan);
        maxLayer = 1;

[~, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
            'omittedWindows', p.omittedWindows, 'Folding', p.Folding, ...
            'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing, ...
            'outlSigma', p.outlSigma, ...
            'EWMA', p.EWMA, ...
            'figuresDir', figuresDir, 'chanName', chan0Name); 
        
% ..st layer

    velmap = rawActmap{1};
    velmap_outl = actmap_outl{1};
    imvelocitymap = imActmap{1}(:, 2:tmax);     %  Imputation (used as an input of computations)
                                                    % Note 2:tmax    

if isempty(MDpixelSize_); error('MD.pixelSize_ is required.'); end
if isempty(MDtimeInterval_); error('MD.timeInterval_ is required.'); end
                                                    

disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])

%%  .txt (export comma delimited files)

fname0 = ['Chan', num2str(iChan)];

dlmwrite(fullfile(figuresDir, [fname0, '_velmap_outl.txt']), velmap_outl, 'precision', 8)
dlmwrite(fullfile(figuresDir, [fname0, 'imvelocitymap.txt']), imvelocitymap, 'precision', 8)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plots 

%%  Raw non-smoothActivityMap prot/act maps
%gaussianFilter = fspecial('gaussian', 13, 3);

%smParam = 1

inputmap = velmap;

%filteredmap = smoothActivityMap(velmap, 'SmoothParam', smParam, 'UpSample', 1);
fvelraw = figure('Visible', p.figFlag);  
figtmp = imagesc(inputmap);
title(chan0Title)
colorbar;colormap(jet) 

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;

%%
saveas(fvelraw, fullfile(figuresDir, ['/raw', fname0, 'Map.png']), 'png')


%%  outl non-smoothActivityMap prot/act maps
%gaussianFilter = fspecial('gaussian', 13, 3);

%smParam = 1

inputmap = velmap_outl;

%filteredmap = smoothActivityMap(velmap, 'SmoothParam', smParam, 'UpSample', 1);
fvelraw = figure('Visible', p.figFlag);  
figtmp = imagesc(inputmap);
title(chan0Title)
colorbar;colormap(jet) 

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;

%%
saveas(fvelraw, fullfile(figuresDir, ['outl_', fname0, 'Map.png']), 'png')



%%  Run:
%%  Boxplots after truncation and outlier detection

%
velBoxTime = figure('Visible', p.figFlag);
boxplot(velmap_outl, 'whisker', Inf);
title(chan0Title)
xlabel('time frame')
set(gca, 'XTick', 10:10:tmax)
set(gca, 'XTickLabel', 10:10:tmax)
 

velBoxWin = figure('Visible', p.figFlag);  
boxplot(velmap_outl', 'whisker', Inf);
title(chan0Title)
xlabel('Window')
set(gca, 'XTick', 10:10:wmax)
set(gca, 'XTickLabel', 10:10:wmax)

%%
saveas(velBoxTime, fullfile(figuresDir, [fname0, 'BoxTime.png']), 'png')
saveas(velBoxWin, fullfile(figuresDir, [fname0, 'BoxWin.png']), 'png')
 

%%  Histogram of velocity. see CV= sm.std/sm.mean*100

sm = summary(velmap_outl(:));
title2 = ['m=', sprintf('%0.2f', sm.mean), ' std=', sprintf('%0.2f', sm.std)];
velHist = figure('Visible', p.figFlag);
histogram(velmap_outl(:));
title({chan0Name, title2});

tmp = velmap_outl(:);

vpos = tmp(tmp >=0);
vneg = tmp(tmp < 0);

sm = summary(vpos);
title2 = ['m=', sprintf('%0.2f', sm.mean), ' std=', sprintf('%0.2f', sm.std)];
vposHist = figure('Visible', p.figFlag);  histogram(vpos);
title({['Positive ', chan0Name], title2})

sm = summary(vneg);
title2 = ['m=', sprintf('%0.2f', sm.mean), ' std=', sprintf('%0.2f', sm.std)];
vnegHist = figure('Visible', p.figFlag);  histogram(vneg);
title({['Negative ', chan0Name], title2})

%%
saveas(velHist, fullfile(figuresDir, [fname0, 'Hist.png']), 'png')
saveas(vposHist, fullfile(figuresDir, [fname0, 'PosHist.png']), 'png')
saveas(vnegHist, fullfile(figuresDir, [fname0, 'NegHist.png']), 'png')


%% topomap topographMD
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topovec = mean(velmap_outl, 2, 'omitnan');
topoMap = nan(wmax, nBandMax_);
topoMap(:, 1) = topovec;

title0 = chan0Title;

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

%%
saveas(topomapFig, fullfile(figuresDir, [fname0, 'topomapFig.png']), 'png')

end

%%  smoothActivityMap prot/act maps

%smParam = 0.8;

inputmap = velmap_outl;

filteredmap = smoothActivityMap(inputmap, 'SmoothParam', p.smParam, 'UpSample', 1);
fvel = figure('Visible', p.figFlag);
figtmp = imagesc(filteredmap);
title(chan0Title)
colorbar;colormap(jet)

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;


%% 
saveas(fvel, fullfile(figuresDir, [fname0, 'Map.png']), 'png')

 
%%  Means plot

inputmap = velmap_outl;

% along time frame
cmeans = mean(inputmap, 1, 'omitnan');
cmeansZ = nanZscore(cmeans) ;

timeFr = 1:tmax;
timeAxis = (timeFr-1)*MDtimeInterval_;

meansTime = figure('Visible', p.figFlag);
plot(timeAxis, cmeansZ);
xlabel('Time (s)'); ylabel('Z-score')
title('Means plot')
legend(chan0Name, 'Location','northoutside','Orientation','horizontal')

% along windows
cmeans = mean(inputmap, 2, 'omitnan');
cmeansZ = nanZscore(cmeans) ;

winIndex = 1:wmax;

meansWin = figure('Visible', p.figFlag);
plot(winIndex, cmeansZ);
title('Means plot')
xlabel('Window'); ylabel('Z-score')
legend(chan0Name, 'Location','northoutside','Orientation','horizontal')


%% 
saveas(meansTime, fullfile(figuresDir, ['/meansTime', fname0, '.png']), 'png')
saveas(meansWin, fullfile(figuresDir, ['/meansWin', fname0, '.png']), 'png')



%%  Coefficient of variation (sd/|mean|)

map = imvelocitymap;

% temporal
stdTemp = std(map, [], 2, 'omitnan');
%meanVec = mean(map, 2, 'omitnan');
%temporalCV = stdVec./abs(meanVec)

% spatial
stdSpat = std(map, [], 1, 'omitnan');

stdFull = [stdTemp ; stdSpat'];

tmax_ = size(map, 2); 
wmax_ = size(map, 1);

grChar = cell(wmax_+tmax_, 1);
for k = 1:wmax_ 
    grChar{k} = 'Temporal';
end
for k = 1:tmax_
    grChar{wmax_+k} = 'Spatial';
end

fvariation = figure('Visible', p.figFlag);
boxplot(stdFull, grChar)
title('Standard deviation (nm/sec)')

%%
saveas(fvariation, fullfile(figuresDir, ['CV', fname0, '.png']), 'png')



%%  temporal SD plot

SDlayers = cell(1, maxLayer);
CVlayers = cell(1, maxLayer);

for indL = 1:maxLayer

    map = actmap_outl{indL};
    % 
    stdTemp = std(map, [], 2, 'omitnan');
    meanTemp = mean(map, 2, 'omitnan');
    temporalCV = stdTemp./abs(meanTemp);

    SDlayers{indL} = stdTemp;
    CVlayers{indL} = temporalCV;
end


%%

fSDlayers = figure('Visible', p.figFlag);

plot(cell2mat(SDlayers))
xlabel('Window'); 
title('Standard deviation')
%legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');
a = 1:maxLayer; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'L'});
    legend(d, 'Location','northoutside','Orientation','horizontal');

ax = gca;
    ax.FontSize = 15;
    

fCVlayers = figure('Visible', p.figFlag);

plot(cell2mat(CVlayers))
xlabel('Window'); 
title('Coeff. Variation')
%legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');
a = 1:maxLayer; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'L'});
    legend(d, 'Location','northoutside','Orientation','horizontal');

ax = gca;
    ax.FontSize = 15;
    

%% save
saveas(fSDlayers, fullfile(figuresDir, ['SDlayers_', fname0, '.png']), 'png')
saveas(fCVlayers, fullfile(figuresDir, ['CVlayers_', fname0, '.png']), 'png')



%% topomap topographMD
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    %topoMap(:, indL) = mean(actmap_outl{indL}, 2, 'omitnan');
    topoMap(:, indL) = SDlayers{indL};
end

title0 = [chan0Title, '-SD'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topographSD_', fname0, '.png']), 'png')





%% 2017/11/04 update

histMeans = figure('Visible', p.figFlag);
histogram(topoMap(:), 'BinMethod', 'fd');
title([chan0Title, '-SDs'])

    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/histSDs_', fname0, '.png']), 'png')

end


%%
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    %topoMap(:, indL) = mean(actmap_outl{indL}, 2, 'omitnan');
    topoMap(:, indL) = CVlayers{indL};
end

title0 = [chan0Title, '-CV'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topographCV_', fname0, '.png']), 'png')




%% 2017/11/04 update

histMeans = figure('Visible', p.figFlag);
histogram(topoMap(:), 'BinMethod', 'fd');
title([chan0Title, '-CVs'])

    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/histCVs_', fname0, '.png']), 'png')

end

%% L2normlayers

%% acmap & LB test, 2017/12/05

map = imvelocitymap;

%%
 
xx1 = 1:ceil(size(map, 2)/2);

%%%% over time
acmapThresh = 2/sqrt(size(map, 2));
acmap = autoCorrMap(map, 'maxLag', xx1(end));  
 
%
%corAvg_autocor = mean(acmap(:, xx1), 1, 'omitnan');
%
acmap1 = acmap(:, 2:size(acmap, 2));
%%%% subplots

acmapOf_ = figure('Visible', p.figFlag);

figtmp = imagesc(acmap1, [-1 1]);
colorbar;colormap(jet)
figtmp.AlphaData = 1-isnan(acmap1);

axis xy;xlabel('Time lag (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick;
ax.XTickLabel = curTick*MDtimeInterval_;

set(gca, 'XGrid', 'on')
title(['AutoCorr of ', chan0Name, ' (Thresh: ', num2str(acmapThresh, 4), ')'])

ax.FontSize = 15;


%%
LBpvec = nan(wmax, 1);
tmp = nan(wmax, 1);
for w = 1:wmax
    if all(~isnan(imvelocitymap(w,:)))
        [~,LBpvec(w), tmp(w)] = lbqtest(imvelocitymap(w, :));
    end
end

LBsig = -log10(LBpvec);
LBsig(~isfinite(LBsig)) = max([LBsig(isfinite(LBsig)); 10]);

fLB = figure('Visible', p.figFlag);
plot(LBsig, 'o-')
refline([0, -log10(0.05)])
xlabel('Window'); ylabel('-log10(P-value)')
LBsig(LBsig > -log10(0.05)) = nan;

hold on
p2 = plot(LBsig, 'o-', 'Color', 'r');

title(['LB test: Significance of autocorr of ', chan0Name])
legend(p2, 'White-noise (Quiescent windows)')

ax = gca;
ax.FontSize = 15;

%%

saveas(acmapOf_, fullfile(figuresDir, ['acmap_', fname0, '.png']), 'png')
saveas(acmapOf_, fullfile(figuresDir, ['acmap_', fname0, '.fig']), 'fig')

saveas(fLB, fullfile(figuresDir, ['LB_significance_', fname0, '.png']), 'png')


%%  ind of LBpvec, threshold=0.05

smLBpvec = movmean(LBpvec, p.movmeanNum);
indActive = (smLBpvec < 0.05);



%%
if strcmp(p.topograph, 'on')

    LBsig2 =  -log10(LBpvec);
    LBsig2(~isfinite(LBsig2)) = max([LBsig2(isfinite(LBsig2)); 10]);
    
iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    %topoMap(:, indL) = mean(actmap_outl{indL}, 2, 'omitnan');
    topoMap(:, indL) = LBsig2;
end

title0 = [chan0Name, ' - LB significance (-log10(P-value))'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['topograph_LB', fname0, '.png']), 'png')
saveas(topomapFig, fullfile(figuresDir, ['topograph_LB', fname0, '.fig']), 'fig')

end


%%
if strcmp(p.topograph, 'on')

    LBsig3 = LBsig2;
    LBsig3(~indActive) = nan;
    
iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    %topoMap(:, indL) = mean(actmap_outl{indL}, 2, 'omitnan');
    topoMap(:, indL) = LBsig3;
end

title0 = [chan0Name, ' - LB significance (-log10(P-value))'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['topograph_LBactive', fname0, '.png']), 'png')
saveas(topomapFig, fullfile(figuresDir, ['topograph_LBactive', fname0, '.fig']), 'fig')

end




%%  Quiescent windows=transparent. smoothActivityMap prot/act maps

%smParam = 0.9;

inputmap = imvelocitymap;

filteredmap = smoothActivityMap(inputmap, 'SmoothParam', p.smParam, 'UpSample', 1);
fvel2 = figure('Visible', p.figFlag);
figtmp = imagesc(filteredmap);
title([chan0Name, ' imputed, classified'])
colorbar;colormap(jet)

alphamap0 = ones(size(inputmap));
alphamap0(~indActive, :) = 0.4;
alphamap0(isnan(inputmap)) = 0; 
figtmp.AlphaData = alphamap0;
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;
    ax.FontSize = 15;
    

%% 
saveas(fvel2, fullfile(figuresDir, [fname0, 'Map_Active.png']), 'png')
saveas(fvel2, fullfile(figuresDir, [fname0, 'Map_Active.fig']), 'fig')

%%  save index

save(fullfile(figuresDir, 'indActive_windowIndex.mat'), 'indActive')


%% avgACF

corAvg_full = mean(acmap(:, xx1), 1, 'omitnan');
corAvg_active = mean(acmap(indActive, xx1), 1, 'omitnan');

favgacf = figure('Visible', p.figFlag);
tlag = (xx1-1)*MDtimeInterval_;
p1 = plot(tlag, corAvg_full);

hold on
p2 = plot(tlag, corAvg_active, '--o');

xlabel('Time lag (s)')
ylabel('Avg. auto correlation')

title('Auto-correlation function')
h = refline(0, 0);
h.Color = [.5 .5 .5];

legend([p1, p2], 'All windows', 'Active windows'); 
%lgd.AutoUpdate = 'off';

ax=gca; ax.FontSize = 15;


%%
saveas(favgacf, fullfile(figuresDir, [fname0, '_avgACF.png']), 'png')
saveas(favgacf, fullfile(figuresDir, [fname0, '_avgACF.fig']), 'fig')

%% output
save(fullfile(figuresDir, 'Vel_meanACF_corAvg_active.mat'), 'corAvg_active', 'corAvg_full', ...
                                'MDtimeInterval_')


%%  checkWindowJump, 2019/03/25

winTrack = checkWindowJump(MD, 'figFlag', p.figFlag);

%
saveas(winTrack, fullfile(figuresDir, 'window1Trajectory.png'), 'png');


%%
disp('====End of mapDescriptives_Vel_LB====')
    
end


