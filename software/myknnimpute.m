function imputed = myknnimpute(mat)
% myknnimpute Impute missing values in a matrix by using the most similar column
% vector (knnimpute.m). If all rows contain NaN, adjust to it.
%
% Jungsik Noh, 2016/01
% Updated, J Noh, 2017/08/28
% J Noh, 2020/01/04. If knnimpute makes an error, it returns all nan.
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

% non-NaN proportion
cind = (sum(~isnan(mat)) ./ size(mat, 1)) > 0.95;

if sum(cind) > 0
    mat1 = mat(:, cind);
    
    % knnimpute.m
    try
        immat1 = knnimpute(mat1);
    catch
        warning('knnimpute() returned an error, so it returns all NaNs.')
        immat1 = nan(size(mat1));
    end    
    imputed = nan(size(mat));
    imputed(:, cind) = immat1;
else
    disp('myknnimpute failed and returned all NaNs.')
    imputed = nan(size(mat));
end

%imputed(:, ~cind) = mat(:, ~cind);

end
