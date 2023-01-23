function permXcorrMean = timePermDistXcorrMean(Actmap1, Actmap2, numPerm, lagMax)
% timePermDistXcorrMean Compute time-wise permuation distribution of a xcorr
% curve using parfor by using nanXcorrMaps.m function.
%
% Jungsik Noh, 2016/10/21
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

% permute time points
tlength = size(Actmap1, 2);
permXcorrMean = nan(numPerm, 2*lagMax+1);
%
parfor i = 1:numPerm;
    permInd = randsample(tlength, tlength);
    
    chan1map = Actmap1(:, permInd);
    chan2map = Actmap2;
        
    xcorrMat_i  = nanXcorrMaps(chan1map, chan2map, 'lagMax', lagMax);
    xcorrMean_i = nanmean(xcorrMat_i, 1);
    
    permXcorrMean(i, :) = xcorrMean_i;
end

