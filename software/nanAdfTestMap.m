function [pvec, adfMap, indicStationairy] = nanAdfTestMap(map, mapName, smParam)
% nanAdfTestMap Perform an augmented Dickey-Fuller (ADF) test to check the
% stationarity of time series of each window and Plot the result.
%
% Output:
%       pvec        - P-values. The null hypothesis is that the time series
%                   is non-stationary. Thus, a small P-value indicates
%                   stationarity.
%       adfMap      - a matlab figure object in which nonstationary windows
%                   are shown with no transparency, and stationary windows
%                   are shown with transparency. 
%
% Updates:
% J Noh, 2018/11/14. Add new output, the indicator of Stationarity.
% J Noh, 2017/10/13. The model option in adftest changed to 'ARD' instead
% of using 'TS' (trend-stationarity) to declare 'non-stationary' when there
% is trend.
% Jungsik Noh, 2016/10/05
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


figFlag = 'off';

indNanRows = any(isnan(map'))';
map2 = map(~indNanRows,:);

%map2 = knnimpute(map')';

% adf
tmax = size(map2, 2);
% automatic lag selection
lag0 = floor((tmax-1)^(1/3));
pvec = NaN*ones(size(map2, 1), 1);

for w = 1:size(map2, 1)
    x = map2(w, :);
    %[h, pval, stat] = adftest(x, 'model', 'TS', 'lags', lag0);
    [~, pval, ~] = adftest(x, 'model', 'ARD', 'lags', lag0);
    %[~, pval, ~] = adftest(x, 'lags', lag0);
    pvec(w) = pval;
end

adfIndex = (pvec < 0.1);           % stationary: 1, nonstationary: 0

% Transform to original map's size
adfIndexOrig = zeros(size(map, 1), 1);
adfIndexOrig(~indNanRows) = adfIndex;
indicStationairy = nan(size(map, 1), 1);
indicStationairy(~indNanRows) = adfIndex;

alphaIndex = (1 - 0.8*adfIndexOrig);    % Nonstationary -> alphaIndex=1 (shown)
%
alphaIndex(indNanRows) = 0;

filteredmap = smoothActivityMap(map, 'SmoothParam', smParam, 'UpSample', 1);
adfMap = figure('Visible', figFlag);
figtmp = imagesc(filteredmap);
colorbar;colormap(jet)
figtmp.AlphaData = repmat(alphaIndex, 1,  tmax);
axis xy;xlabel('Time frame');ylabel('Window')
title({mapName, ['total num: ', num2str(size(map2,1)), ' num of nonstationary TS: ', num2str(sum(~adfIndex))] } ) 

end


% ftn:end
