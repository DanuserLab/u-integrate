function sigma = long_run_variance(data,band)
% Returns the long run variance of data, a T*m matrix.
%
% INPUTS
%    data      [double]  T*m matrix, where T is the number of observations and m the number of variables.
%    band      [double]  scalar, the bandwidth parameter.
%
% OUTPUTS
%    sigma     [double]  m*m matrix.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2009-2010 Dynare Team
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
 
 verbose = 1;
 
 if nargin<2
     [T,m] = size(data);
     band = ceil(T^(1/4));
     if verbose
         disp(['Bandwidth parameter is equal to ' num2str(band) '.'])
     end
 end
 
 gamma = multivariate_sample_autocovariance(data,band);
 sigma = gamma(:,:,1);
 for i=1:band
     sigma = sigma + bartlett(i,band)*(gamma(:,:,i+1)+transpose(gamma(:,:,i+1)));
 end
 
 function w = bartlett(i,n)
 w = 1 - i / (n+1);

