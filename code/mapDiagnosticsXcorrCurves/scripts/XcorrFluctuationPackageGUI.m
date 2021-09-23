function varargout = XcorrFluctuationPackageGUI(varargin)
% Launch the GUI for the XcorrFluctuation Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Qiongjing (Jenny) Zou, Sep 2018
%

  % Comment out below, b/c I do not want to switch to the MDs as packageGUI's input, when use ML as input.
% if nargin>0 && isa(varargin{1},'MovieList')
%     varargout{1} = packageGUI('XcorrFluctuationPackage',[varargin{1}.getMovies{:}],...
%         varargin{2:end}, 'ML', varargin{1});
% else
    varargout{1} = packageGUI('XcorrFluctuationPackage',varargin{:});
% end

end