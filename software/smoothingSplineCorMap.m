function fittedCorVec = smoothingSplineCorMap(corMap)
% smoothingSplineCorMap Compute an averaged correlation curve over windows
% by applying smoothingspline via fit.m function. The smoothingParam is
% chosen to be twice of the default one. This function handles NaN
% correlations by simply removing in the fit.
% Jungsik Noh, 2017/03/31
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

xl = 1:size(corMap, 2);
yy = corMap(:);
tmpA = repmat(xl, size(corMap, 1), 1);
xx = tmpA(:);

% remove nan's
xx = xx(~isnan(yy));
yy = yy(~isnan(yy));

[~,~,tmpfit] = fit(xx, yy, 'smoothingspline');
smParSp = min(1, tmpfit.p*2);   % make the fit less smooth than the optimal

fitsp = fit(xx, yy, 'smoothingspline', 'SmoothingParam', smParSp);
fittedCorVec = fitsp(xl);

fittedCorVec = reshape(fittedCorVec, 1, []);

end
