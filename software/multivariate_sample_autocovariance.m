function autocov = multivariate_sample_autocovariance(data,q)
% Computes the autocovariance function of multivariate time series.
%
%
% INPUTS
%   data            [double]       T*m matrix of data.
%   q               [integer]      Order of the autocovariance function.
%
% OUTPUTS
%   autocov         [double]       m*m*(q+1) array, autocovariance function.
%
% SPECIAL REQUIREMENTS
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

% Copyright (C) 2009-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

[T,m] = size(data);

autocov = zeros(m,m,q+1);

demeaned_data = data - repmat(mean(data),T,1);

for i = 0:q
    autocov(:,:,i+1) = transpose(demeaned_data(1:T-i,:))*demeaned_data(i+1:T,:);
end

autocov = autocov/T;