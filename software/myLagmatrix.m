function laggedTScolMat = myLagmatrix(TScolMat, lagOrder)
% lagged TS column matrix with intial NaN's imputed by means
%
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


X = TScolMat;

%%
Xlmat = lagmatrix(X, 1:lagOrder);
Xlmat1 = Xlmat;

% nan impute with mean
for l = 1:lagOrder
    for j = 1:l
        Xlmat1(j, 1+(l-1)*size(X, 2):l*size(X,2)) = mean(X, 1, 'omitnan');
    end
end

laggedTScolMat = Xlmat1;

end
