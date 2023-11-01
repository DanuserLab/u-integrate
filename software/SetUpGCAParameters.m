function SetUpGCAParameters(movieListOrProcess, varargin)
% SetUpGCAParameters specifies subsequent parameters from the first GCA parameters. 
%
% Jungsik Noh, 7/2021
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
 
%% Input
ip = inputParser;
ip.addRequired('movieListOrProcess', @isProcessOrMovieList);
ip.addOptional('param',[], @isstruct);
ip.parse(movieListOrProcess, varargin{:});
p = ip.Results;
paramsIn = p.param;

%% set up subsequent parameters

[movieList, process] = getOwnerAndProcess(movieListOrProcess,'SetUpGCAParametersProcessML',true);
p = parseProcessParams(process, paramsIn);  % If parameters are explicitly given, they should be used
                                        % rather than the one stored in XcorrAnalysisProcessML

% p.avgLengthOfProtRetCycle is required for setting up multiple pars like 
% movMedFrameSize, ARorder, lagMax, etc                                         
if isempty(p.avgLengthOfProtRetCycle) || ~isfinite(p.avgLengthOfProtRetCycle)
    error('avgLengthOfProtRetCycle (an integer for frames) is required.')
end
                                        
% calculate internal parameters from the global parameters like
% chanCodeWithPreprocess = 21 (lowFrequency Subtracted chan1)
% Low-Frequency-Subtracted signals, chanCodeWithPreprocess1 = 2x
if p.lowFreqSubtraction
    p.chanCodeWithPreprocess = 20 + p.ChannelIndex;
    p.chanNameWithPreprocess = ...
        {['LFS.', p.ChannelName{1}], ['LFS.', p.ChannelName{2}]};
else
    p.chanCodeWithPreprocess = p.ChannelIndex;
    p.chanNameWithPreprocess = p.ChannelName;
end
            
    % For other various preprocessing methods in
    % mapOutlierImputation.m, directly specify the
    % chanCodeWithPreprocess1 and chanNameWithPreprocess1.
    % e.g., low-Frequency of signals, chanCodeWithPreprocess1 = 3x

% Update parameters

SetUpGCAParametersProcessMLid = movieList.getProcessIndex('SetUpGCAParametersProcessML');

parseProcessParams(movieList.processes_{SetUpGCAParametersProcessMLid}, p);

end
