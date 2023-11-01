function MLsummary_quiescentWindow(ML, LBdirName, varargin)
% MLsummary_quiescentWindow
%
% Output: .png/.fig/.mat files are saved in the 'ML.outputDirectory_/MLdiagnosticPlots' by default. 
%
% Updated:
% J Noh, 2019/06/11. Add output of active_Chan0 maps and acmap_Chan0.
% J Noh, 2019/03/25. checkWindowJump added.
% J Noh, 2019/02/27. Compute ACF summary and lag (s) of minimum meanACF.
% Updated by Qiongjing (Jenny) Zou Nov 2018
% Added new parameter called 'outDirName'
%
% J Noh, 2018/11/13
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

%% parameter 

ip = inputParser;
ip.addParameter('ratioThreshold', 0.5);
ip.addParameter('outDirName', 'MLdiagnosticPlots');

parse(ip, varargin{:});
p = ip.Results;

% outDir
outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%% load ML, example MD
ML.getMovies()
md1 = ML.getMovie(1);
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])

%%
MDs = ML.getMovies();
num = numel(MDs);
indActive2 = cell(num, 1);
cellLabels = cell(num, 1);
acfvecsize = nan(num, 1);
acfArr = cell(num, 1);

for i=1:num
    fname = fullfile(MDs{i}.outputDirectory_, LBdirName, 'indActive_windowIndex.mat');
    S = load(fname);
    indActive2{i} = S.indActive;
    %
    [~, cellLab0, ~] = fileparts(MDs{i}.outputDirectory_);
    cellLabels{i} = [cellLab0(1:end)]; 
    % for ACF summary
    S2 = load(fullfile(MDs{i}.outputDirectory_, LBdirName, 'Vel_meanACF_corAvg_active.mat')); % Avg_autocor
    acfvecsize(i) = size(S2.corAvg_active, 2);
    acfArr{i} = S2.corAvg_active;
end
% MDtimeInterval_ from the last MD
MDtimeInterval_ = S2.MDtimeInterval_;

ratioVec = nan(1, num);
for i=1:num
    ratioVec(i) = 1 - nansum(indActive2{i})/sum(~isnan(indActive2{i}));
end

%%
fLBbar = figure;
bar(1:num, ratioVec)

set(gca, 'XTick', 1:num)
set(gca, 'XTickLabel', cellLabels)
set(gca, 'XTickLabelRotation', 45)
xlabel('movieData index')
ylabel('Proportion of quiescent windows')
ylim([0,1])
h = refline([0, p.ratioThreshold]); h.Color = 'r';
title({'Remove movieData for which quiescent windows are '; ['more than ', num2str(round(100*p.ratioThreshold)), '%']})

%% saveas
saveas(fLBbar, fullfile(outDir, 'percentageOfQuiescentWindowsPerMD.png'), 'png')
saveas(fLBbar, fullfile(outDir, 'percentageOfQuiescentWindowsPerMD.fig'), 'fig')
save(fullfile(outDir, 'ratioOfQuiescentWindowsPerMD.mat'), 'ratioVec')

%% acf summary

acfvecsize0 = min(acfvecsize);
acfArr2 = nan(num, acfvecsize0);
for i=1:num
    acfArr2(i,:) = acfArr{i}(1:acfvecsize0);
end
MLmeanACF = mean(acfArr2, 1, 'omitnan');
[~, minid] = min(MLmeanACF);
%
facfsummary = figure;
tLag = [0:acfvecsize0-1].* MDtimeInterval_;
acfminTime = tLag(minid);
    
    p1 = plot(tLag, MLmeanACF, 'r');
    p1.LineWidth = 2;    
    hold on
    
    plot(tLag, acfArr2, 'Color', 'k', 'LineWidth', 0.5)
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    
    title(['Vel (ACF): meanACF min at ', num2str(acfminTime), ' s'] )
    
    xlabel('Time lag (s)');ylabel('Correlation')
    set(gca, 'XGrid', 'on')
    ax = gca;
    ax.FontSize = 14;
    
%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
%ylim([-1, 1])
xlim([min(tLag), max(tLag)])

%% saveas

saveas(facfsummary, fullfile(outDir, ['Vel_LB', '_acf', '.png']), 'png')
saveas(facfsummary, fullfile(outDir, ['Vel_LB', '_acf', '.fig']), 'fig')

save(fullfile(outDir, 'Vel_LB_acfs.mat'), 'acfArr2', 'MLmeanACF', 'acfminTime', 'MDtimeInterval_')


%%  summarize checkWindowJump

for i=1:num
    fname2 = fullfile(MDs{i}.outputDirectory_, LBdirName, 'window1Trajectory.png');
    copyfile(fname2, fullfile(outDir, [cellLabels{i}, '_window1Trajectory.png']) );    
end

%% 
for i=1:num
    fname3 = fullfile(MDs{i}.outputDirectory_, LBdirName, 'Chan0Map_Active.png');
    copyfile(fname3, fullfile(outDir, [cellLabels{i}, '_Chan0Map_Active.png']) );    
end

%% acmap_Chan0
for i=1:num
    fname4 = fullfile(MDs{i}.outputDirectory_, LBdirName, 'acmap_Chan0.png');
    copyfile(fname4, fullfile(outDir, [cellLabels{i}, 'acmap_Chan0.png']) );    
end



end
