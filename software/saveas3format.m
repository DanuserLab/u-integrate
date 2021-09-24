function saveas3format(fig, outputDir, fname)
% saveas3format SAVE a figure in .fig, .eps, and .pdf formats.
%
% Jungsik Noh, 2018/01/23
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

fullname1 = fullfile(outputDir, [fname, '.fig']);
fullname2 = fullfile(outputDir, [fname, '.eps']);
fullname3 = fullfile(outputDir, [fname, '.pdf']);

saveas(fig, fullname1, 'fig');
saveas(fig, fullname2, 'epsc');
saveas(fig, fullname3, 'pdf');


end
