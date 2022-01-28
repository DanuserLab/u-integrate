function acfMat1 = autoCorrMap(map, varargin)
% autoCorrMap Compute auto-correlations over time for each window (over
% columns for each row) while handling NaN's.
% Jungsik Noh, 2016/10/29
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

tmax = size(map, 2);
wmax = size(map, 1);
maxLag0 = floor(tmax/4);

ip = inputParser;
ip.addRequired('map', @ismatrix)
ip.addParameter('maxLag', maxLag0);

ip.parse(map, varargin{:});
maxLag = ip.Results.maxLag;

acfMat1 = zeros(wmax, maxLag+1);

for w = 1:wmax
    acf0 = nanXcorrelation(map(w, :), map(w, :), maxLag);   % nanXcorrelation's ftn is used.
    acf = acf0(maxLag+1:end);
    acfMat1(w, :) = acf;
end

