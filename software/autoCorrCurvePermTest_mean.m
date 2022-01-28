function [acCurve, Avg_autocor] = autoCorrCurvePermTest_mean(map, mapName, MDtimeInterval_, ...
                    numPerm, parpoolNum, rseed, varargin)
% autoCorrCurvePermTest Draw a curve of the temporal auto correlation means 
% together with its confidence bound under the null.
% Jungsik Noh, 2016/10/19
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

% Updated 2017/03/31: 
%           Instead of mean(), it uses smoothingSpline.

ip = inputParser;
ip.addParameter('figFlag', 'off');
parse(ip, varargin{:})
p = ip.Results;
figFlag = p.figFlag;

%% due to parfor
if isempty(gcp('nocreate')); parpool('local', parpoolNum); end
% rng set up
rng('default'); rng(rseed)


%%  auto corr map/curve

wmax = size(map, 1);
maxLag = ceil(size(map, 2)/2); 

xx1 = 1:maxLag;
acmapThresh = 2/sqrt(size(map, 2));

% centered map (minus rowmeans)
map = map - mean(map, 2, 'omitnan')*ones(1, size(map,2));

acmap = autoCorrMap(map, 'maxLag', xx1(end));  
 
% if using mean:
Avg_autocor = mean(acmap, 1, 'omitnan');
% if using smoothingspline
%if all(isnan(acmap(:)))
%    Avg_autocor = nan(1, size(acmap, 2));
%else
%    Avg_autocor = smoothingSplineCorMap(acmap);
%end



%% time-wise permutated corr means

tlength = size(map, 2);
permAcorrMeans = nan(numPerm, maxLag+1);
%
permInd = cell(numPerm, 1);
permMap = cell(numPerm, 1);
for i = 1:numPerm;
    permInd{i} = randsample(tlength, tlength);
    permMap{i} = map(:, permInd{i});
end

parfor i = 1:numPerm;
%    acmap_i  = autoCorrMap(permMap{i}, 'maxLag', maxLag);
    
    acmap_i = zeros(wmax, maxLag+1);
    for w = 1:wmax
        xcr = xcorr(map(w, :), permMap{i}(w, :), maxLag, 'coeff');   % ... is used.
        xcr(1:maxLag) = [];
        acmap_i(w, :) = xcr;
    end
    
    % autocor mean curve per permutation
    acmean_i = mean(acmap_i, 1, 'omitnan');
    
    permAcorrMeans(i, :) = acmean_i;
end

numRowAcmap = sum(~isnan(sum(acmap, 2)));

%
uCL = quantile(permAcorrMeans, 0.975);
lCL = quantile(permAcorrMeans, 0.025);



%%  curve
%%  autocorr curves with permutation confidence limits under the null

acCurveThresh = acmapThresh/sqrt(numRowAcmap);

tlag = [0, xx1]*MDtimeInterval_;

acCurve = figure('Visible', figFlag);                                  % Figure name

p1 = plot(tlag, Avg_autocor);
xlabel('Time lag (s)')
ylabel('Avg. auto correlation')
title(mapName)
h = refline(0, 0);
h.Color = [.5 .5 .5];
set(gca, 'XGrid', 'on') 
    
hold on  
tLag1 = xx1*MDtimeInterval_;
p2 = plot(tLag1, uCL(2:(maxLag+1)), 'b--');
plot(tLag1, lCL(2:(maxLag+1)), 'b--')

h = refline([0, acCurveThresh]);
h.Color = [.5 .5 .5];
h.LineStyle = '--';
h = refline([0, -acCurveThresh]);
h.Color = [.5 .5 .5];
h.LineStyle = '--';

legend([p1, p2], 'Avg. Auto corr.', 'Confidence Bound')


end