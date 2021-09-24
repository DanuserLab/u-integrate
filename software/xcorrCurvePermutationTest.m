function xcorrMat = xcorrCurvePermutationTest(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
        fName0, MDtimeInterval_, figuresDir, varargin)
% xcorrCurvePermutationTest Draw 
%   (1) cross correlation maps between two channels for each layer, 
%   (2) cross correlation curves, and 
%   (3) their permutation confidence limit of the curve under the null.
% It computes cross correlations in a fashion that can handle many NaN's 
% by utilizing nanXcorrMaps.m function.
%
% Jungsik Noh, 2016/10/21
%
% Copyright (C) 2021, Danuser Lab - UTSouthwestern 
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


% Updated 
% J Noh, 2021/07/21. Due to Java memory error in GUI, add figure close
% commands.
% J Noh, 2019/07/02. colOrd1 is revised.
% J Noh, 2018/11/14. Edit legends for layers.
% J Noh, 2018/11/14. Remove parpoolNum, parpool().
% J Noh, 2018/04/03. Fix remaining nan bugs.
% Updated:
% J Noh, 2018/01/29. Correlation map and curve plots are now saved also in
% .fig format. Legend option is adjusted. 
% 2017/03/31: 
%           Instead of mean(), it uses smoothingSpline.


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
if ~isfolder(figuresDir); mkdir(figuresDir); end

%% due to parfor
%if isempty(gcp('nocreate')); parpool('local', p.parpoolNum); end
% rng set up
rng('default'); rng(p.rseed)


%% xcorrMapPlot
[~, xcmap_fig] =  xcorrMapPlot(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
    MDtimeInterval_, 'smParam', 1, 'lagMax', p.lagMax,  ...
     'fullRange', 1, 'figFlag', p.figFlag);

%% saveas
    for indL = 1:layerMax
        saveas(xcmap_fig{indL}, fullfile(figuresDir, ['xcmap_', fName0, '_', num2str(indL), 'L.png']), 'png')
        pause(0.1)
    end

    
%% xcorrMapPlot with range
[xcorrMat, xcmap_fig] =  xcorrMapPlot(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName,  ...
    MDtimeInterval_, 'smParam', 1, 'lagMax', p.lagMax, ...
    'fullRange', p.fullRange, 'figFlag', p.figFlag);


%% saveas
for indL = 1:layerMax
    saveas(xcmap_fig{indL}, fullfile(figuresDir, ['xcmap_', fName0, '_range_', num2str(indL), 'L.png']), 'png')
    saveas(xcmap_fig{indL}, fullfile(figuresDir, ['xcmap_', fName0, '_range_', num2str(indL), 'L.fig']), 'fig')
    pause(0.1)
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(xcmap_fig{indL})
end

%% time-wise permutation 

permXcorrMean = cell(layerMax, 1);
numRowXcorrMat = cell(layerMax, 1);
uCL = cell(layerMax, 1);
lCL = cell(layerMax, 1);

for indL = 1:layerMax

%    Actmap1 = squeeze(ch1Actmap(:, indL, :));
%    Actmap2 = squeeze(ch2Actmap(:, indL, :));
    Actmap1 = ch1Actmap{indL};
    Actmap2 = ch2Actmap{indL};


    permXcorrMean{indL} = timePermDistXcorrMean(Actmap1, Actmap2, p.numPerm, p.lagMax);  % permute Actmap1

    numRowXcorrMat{indL} = sum(~isnan(sum(xcorrMat{indL}, 2)));

    %

    uCL{indL} = quantile(permXcorrMean{indL}, 0.975);
    lCL{indL} = quantile(permXcorrMean{indL}, 0.025);

end

%%  xcorr curves with permutation confidence limits under the null

%lagGrid = round(p.lagMax/2);
%xcmapXtick = 1:lagGrid:(p.lagMax*2+1);
%xcmapXticklabel = round((-p.lagMax:lagGrid:p.lagMax) * MDtimeInterval_, 2);

lagMax = p.lagMax;
lagGrid = floor(p.lagMax/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeInterval_, 2);

colOrd = get(gca,'ColorOrder');
%colOrd1 = colOrd(1:5, :);          % Revise later together with legend.
colOrd1 = repmat(colOrd, 10, 1);    % Revised setting maxLayer=70.
 
xcorrMean = nan(p.lagMax*2+1, layerMax);
for indL = 1:layerMax
    % if using mean
    % xcorrMean(:, indL) = nanmean(xcorrMat{indL}, 1);
    % if using smoothingspline
    if ~all(isnan(xcorrMat{indL}))
        xcorrMean(:, indL) = smoothingSplineCorMap(xcorrMat{indL});
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
uCLs = cell2mat(uCL);
lCLs = cell2mat(lCL);
 
for indL = 1:layerMax
    plot(uCLs(indL,:)', '--', 'Color', colOrd1(indL,:))
end

for indL = 1:layerMax
    plot(lCLs(indL,:)', '--', 'Color', colOrd1(indL,:))
end

    %legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');
    a = 1:layerMax; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'L'});
    legend(d, 'Location','northoutside','Orientation','horizontal');
    

%% saveas
saveas(xcCurve_xxx_permT, fullfile(figuresDir, ['/xcCurve_', fName0, '_permT.png']), 'png')
saveas(xcCurve_xxx_permT, fullfile(figuresDir, ['/xcCurve_', fName0, '_permT.fig']), 'fig')

%% close figures
pause(0.5) 
close(xcCurve_xxx_permT) 

%%  output: xcorrMat for pooling

%xcorrMat_tmp = xcorrMat;     %% Name
save(fullfile(figuresDir, ['xcorrMat_', fName0, '.mat']), 'xcorrMat')
save(fullfile(figuresDir, ['permXcorrMean_', fName0, '.mat']), 'permXcorrMean')
save(fullfile(figuresDir, ['numRowXcorrMat_', fName0, '.mat']), 'numRowXcorrMat')


end

