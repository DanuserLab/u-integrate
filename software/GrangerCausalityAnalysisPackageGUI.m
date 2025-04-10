function varargout = GrangerCausalityAnalysisPackageGUI(varargin)
% Launch the GUI for the GrangerCausalityAnalysis Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%  
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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

  % Comment out below, b/c I do not want to switch to the MDs as packageGUI's input, when use ML as input.
% if nargin>0 && isa(varargin{1},'MovieList')
%     varargout{1} = packageGUI('GrangerCausalityAnalysisPackage',[varargin{1}.getMovies{:}],...
%         varargin{2:end}, 'ML', varargin{1});
% else
    varargout{1} = packageGUI('GrangerCausalityAnalysisPackage',varargin{:});
% end

end