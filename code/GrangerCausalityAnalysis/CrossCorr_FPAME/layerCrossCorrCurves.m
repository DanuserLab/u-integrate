function xcorrMatArr = layerCrossCorrCurves(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
        fName0, MDtimeInterval_, figuresDir, varargin)
% layerCrossCorrCurves() Draw 
%   (1) cross correlation maps between two channels for each layer and between layers, 
%   (2) cross correlation curves. 
% It computes cross correlations in a fashion that can handle many NaN's 
% by utilizing nanXcorrMaps.m function.
%
% Updated: 
% J. Noh. 2019/12/12. Change smoothingSplineCorMap into nanmean due to
% lag=0 correlations. 
%
% Jungsik Noh, 2019/03/25



%%  initialization

%wmax = size(ch1Actmap, 1);
tmax = size(ch1Actmap{1}, 2);
layerMax = numel(ch1Actmap);

ip = inputParser;
ip.addParameter('figFlag', 'off');
ip.addParameter('lagMax', round(tmax/4), @isnumeric);
ip.addParameter('fullRange', false);
ip.addParameter('parpoolNum', 4);
ip.addParameter('rseed', 'shuffle');
ip.addParameter('numPerm', 1000);

parse(ip, varargin{:})
p = ip.Results;
%figFlag = p.figFlag;
set(groot,'defaultLegendAutoUpdate','off') 

%
if ~isdir(figuresDir); mkdir(figuresDir); end

%% due to parfor
%if isempty(gcp('nocreate')); parpool('local', p.parpoolNum); end
% rng set up
%rng('default'); rng(p.rseed)


%% xcorrMapPlot
[~, xcmap_figArr] =  layerCrossCorrMapPlot(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
    MDtimeInterval_, 'smParam', 1, 'lagMax', p.lagMax,  ...
     'fullRange', 1, 'figFlag', p.figFlag);

%% saveas
    for indL = 1:layerMax
    for indL2 = 1:layerMax
        saveas(xcmap_figArr{indL, indL2}, fullfile(figuresDir, ...
            ['xcmap_', fName0, '_', num2str(indL), 'L-', num2str(indL2), 'L.png']), 'png')
        pause(0.1)
    end
    end

    
%% xcorrMapPlot with range
[xcorrMatArr, xcmap_figArr] =  layerCrossCorrMapPlot(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName,  ...
    MDtimeInterval_, 'smParam', 1, 'lagMax', p.lagMax, ...
    'fullRange', p.fullRange, 'figFlag', p.figFlag);


%% saveas
for indL = 1:layerMax
for indL2 = 1:layerMax    
    saveas(xcmap_figArr{indL, indL2}, fullfile(figuresDir, ...
        ['xcmap_', fName0, '_range_', num2str(indL), 'L-', num2str(indL2), 'L.png']), 'png')
    saveas(xcmap_figArr{indL, indL2}, fullfile(figuresDir, ...
        ['xcmap_', fName0, '_range_', num2str(indL), 'L-', num2str(indL2), 'L.fig']), 'fig')
    pause(0.1)
end
end

xcorrMat = cell(layerMax, 1);
for indL = 1:layerMax
    xcorrMat{indL} = xcorrMatArr{indL, indL};
end


%% time-wise permutation. Dropped. 2019/03.

%%  xcorr curves with permutation confidence limits under the null

%lagGrid = round(p.lagMax/2);
%xcmapXtick = 1:lagGrid:(p.lagMax*2+1);
%xcmapXticklabel = round((-p.lagMax:lagGrid:p.lagMax) * MDtimeInterval_, 2);

lagMax = p.lagMax;
lagGrid = floor(p.lagMax/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeInterval_, 2);

colOrd = get(gca,'ColorOrder'); 
nRow0 = size(colOrd,1); multipleNum = ceil(layerMax/nRow0);
colOrdExt = repmat(colOrd, multipleNum, 1);
colOrd1 = colOrdExt(1:layerMax, :);                   % Revise later together with legend.
 
xcorrMean = nan(p.lagMax*2+1, layerMax);
for indL = 1:layerMax
    % if using mean
    % xcorrMean(:, indL) = nanmean(xcorrMat{indL}, 1);
    if ~all(isnan(xcorrMatArr{indL, indL}))
    % if using smoothingspline        
        %xcorrMean(:, indL) = smoothingSplineCorMap(xcorrMatArr{indL, indL});
        % if using nanmean
        xcorrMean(:, indL) = nanmean(xcorrMatArr{indL, indL});
    end

end

    xcCurve_xxx_permT = figure('Visible', p.figFlag);   %%  name
    plot(xcorrMean)
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, '_t)'])
    xlabel('Time lag (s)');ylabel('Correlation')

    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    refline([0 0])
    refline([Inf 0])
    set(gca, 'XGrid', 'on')
    
    hold on 
%uCLs = cell2mat(uCL);
%lCLs = cell2mat(lCL);

    %legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');
    a = 1:layerMax; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'L'});
    legend(d, 'Location','northoutside','Orientation','horizontal');
    

%% saveas
saveas(xcCurve_xxx_permT, fullfile(figuresDir, ['/xcCurve_', fName0, '_permT.png']), 'png')
saveas(xcCurve_xxx_permT, fullfile(figuresDir, ['/xcCurve_', fName0, '_permT.fig']), 'fig')

%%  xcorr curves Between layers

xcCurve_betnLayer = cell(1, layerMax);
for indL2 = 1:layerMax

xcorrMeanTmp = nan(p.lagMax*2+1, layerMax);
for indL = 1:layerMax
    %Li = min(indL, indL2); Lj = max(indL, indL2);
    if ~all(isnan(xcorrMatArr{indL, indL2}))
    % if using smoothingspline
        %xcorrMeanTmp(:, indL) = smoothingSplineCorMap(xcorrMatArr{indL, indL2});
        % if using nanmean               
        xcorrMeanTmp(:, indL) = nanmean(xcorrMatArr{indL, indL2});
    end

end

    xcCurve_betnLayer{indL2} = figure('Visible', p.figFlag);   %%  name
    plot(xcorrMeanTmp)
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, '_t)', ...
        ' for ', ch2ActmapName, '-', num2str(indL2), 'L'])
    xlabel('Time lag (s)');ylabel('Correlation')
    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    refline([0 0])
    refline([Inf 0])
    set(gca, 'XGrid', 'on')
    
    hold on 
%uCLs = cell2mat(uCL);
%lCLs = cell2mat(lCL);

    %legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');
    a = 1:layerMax; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'L'});
    legend(d, 'Location','northoutside','Orientation','horizontal');

end

%% saveas xcCurve_betnLayer
for indL2 = 1:layerMax
    saveas(xcCurve_betnLayer{indL2}, fullfile(figuresDir, ['/xcCurve_', fName0, '_to', num2str(indL2), 'L.png']), 'png')
    saveas(xcCurve_betnLayer{indL2}, fullfile(figuresDir, ['/xcCurve_', fName0, '_to', num2str(indL2), 'L.fig']), 'fig')
end


%%  output: xcorrMat for pooling

%xcorrMat_tmp = xcorrMat;     %% Name
save(fullfile(figuresDir, ['xcorrMat_', fName0, '.mat']), 'xcorrMat')
save(fullfile(figuresDir, ['xcorrMatArr_', fName0, '.mat']), 'xcorrMatArr')


%save(fullfile(figuresDir, ['permXcorrMean_', fName0, '.mat']), 'permXcorrMean')
%save(fullfile(figuresDir, ['numRowXcorrMat_', fName0, '.mat']), 'numRowXcorrMat')


end

