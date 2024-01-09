function xCorrMap = nanXcorrMaps(xmap, ymap, varargin)
% nanXcorrMaps Compute cross correlations between two activity maps. The cross
% correlations at lag h are Corr(xmap_{t+h}, ymap_t).
% This function can handle many NaN's due to using nanXcorrelatioon.m
% function.
%
% Jungsik Noh, 2016/10/xx
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

%%
tmax0 = size(xmap, 2);

ip = inputParser;
ip.addRequired('xmap',@(x)(ismatrix(x) && min(size(x)) > 1));
ip.addRequired('ymap',@(x)(ismatrix(x) && min(size(x)) > 1));

ip.addParameter('lagMax', round(tmax0/4), @isnumeric);
ip.parse(xmap, ymap, varargin{:});
p = ip.Results;

%%
map1 = xmap' - repmat(mean(xmap', 1, 'omitnan'), size(xmap', 1), 1);
map2 = ymap' - repmat(mean(ymap', 1, 'omitnan'), size(ymap', 1), 1);

%lagMax = 50;
lagMax = p.lagMax;
numWindows = size(map1, 2);

xCorrMap = zeros(numWindows, 2*lagMax+1);

for w=1:numWindows
    %xCorrMap(w, :) = xcorr(map1(:, w), map2(:, w), lagMax, 'coeff');    
    xCorrMap(w, :) = nanXcorrelation(map1(:, w), map2(:, w), lagMax);
end

 

