function CMLags = xcorrMatToCMLag(xcorrMatCell)
% xcorrMatToCMLag Extract a feature of correlation maxima and their lags
% from the xcorrMat per movie.
%   OUTPUT:     r0, lag0 per layer
% 
% Jungsik Noh, 02/07/2017
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
% J Noh, 2018/04/03. Fix all nan case.
% 2017/03/31: 
%           Instead of mean(), it uses smoothingSpline.

numLayer = numel(xcorrMatCell);
lagMax = (size(xcorrMatCell{1}, 2) - 1)/2;

CMLags = cell(numLayer, 1);
for l = 1:numLayer
    xcmat = xcorrMatCell{l};
    %xcmeanVec = mean(xcmat, 1, 'omitnan');
    % all nan case
    if all(isnan(xcmat))
        xcmeanVec = mean(xcmat, 1, 'omitnan');
    else
    xcmeanVec = smoothingSplineCorMap(xcmat);
    end
    
    [r0, i] = max(abs(xcmeanVec));
    l0 = i - lagMax - 1;
    r1 = xcmeanVec(i);
    CMLags{l} = [r1, l0];
end


end
