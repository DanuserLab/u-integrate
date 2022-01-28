function MLsummary_ADFtest(ML, mapDescDirName, iChan, chanName, varargin)
% MLsummary_ADFtest
%
% Output: .png/.fig/.mat files are saved in the 'ML.outputDirectory_/MLdiagnosticPlots' by default. 
%
% J Noh, 2018/11/13
%
% Copyright (C) 2022, Danuser Lab - UTSouthwestern 
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

% Updated by Qiongjing (Jenny) Zou Nov 2018
% Added new parameter called 'outDirName'


%% parameter 

ip = inputParser;
ip.addParameter('outDirName', 'MLdiagnosticPlots');

parse(ip, varargin{:});
p = ip.Results;

% outDir
outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%% load ML, example MD
ML.getMovies();
md1 = ML.getMovie(1);
%[~, cellLab0, ~] = fileparts(md1.outputDirectory_);
%disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])

fname0 = ['Chan', num2str(iChan)];
fname = fullfile(md1.outputDirectory_, mapDescDirName, ['indicStationarity_', fname0, '.mat']);
S = load(fname);
indTmp = S.indicStationarity;
maxLayer = numel(indTmp);

%%
MDs = ML.getMovies();
num = numel(MDs);
indicStationarity2 = cell(num, 1);
colLabels = cell(num, 1);
ratioMat = nan(num, maxLayer);

for i=1:num
    fname = fullfile(MDs{i}.outputDirectory_, mapDescDirName, ['indicStationarity_', fname0, '.mat']);
    S = load(fname);
    indicStationarity2{i} = S.indicStationarity;
    %
    for indL = 1:maxLayer
        indStat = indicStationarity2{i}{indL};
        rtmp = 1 - nansum(indStat)/sum(~isnan(indStat));
        ratioMat(i, indL) = rtmp;
        colLabels{i} = ['MD', num2str(i)];
    end
end


%%
fADFbar = figure;
if (num == 1)
    bar(ratioMat)
    a = 1:maxLayer; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'layer'}); %e = reshape(d, 1, []);
    set(gca, 'XTickLabel', d)
else
    bar(1:num, ratioMat)
    a = 1:maxLayer; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'layer'}); %e = reshape(d, 1, []);
    legend(d, 'Location', 'northoutside', 'Orientation', 'horizontal')
    set(gca, 'XTick', 1:num)
    set(gca, 'XTickLabel', colLabels)
    
end

xlabel('movieData Index')
ylabel({'Proportion of windows '; 'with systematic trends'})
ylim([0,1])
h = refline([0, 0.5]); h.Color = 'r';
title({'Time series with systematic trends '; 'are known to produce spurious cross correlations!'; chanName})
%title({'Remove movieData for which quiescent windows are '; ['more than ', num2str(round(100*p.ratioThreshold)), '%']})

%a = 1:maxLayer; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'layer'}); %e = reshape(d, 1, []);
%legend(d, 'Location', 'northoutside', 'Orientation', 'horizontal')


%% saveas
saveas(fADFbar, fullfile(outDir, ['percentageOfWindowsWithTrendsPerMD_', fname0, '.png']), 'png')
saveas(fADFbar, fullfile(outDir, ['percentageOfWindowsWithTrendsPerMD_', fname0, '.fig']), 'fig')
save(fullfile(outDir, ['ratioMatOfWindowsWithTrendsPerMD_', fname0, '.mat']), 'ratioMat')

%%

disp('==== MLsummary_ADFtest is finished!! ====')

end
