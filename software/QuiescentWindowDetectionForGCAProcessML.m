classdef QuiescentWindowDetectionForGCAProcessML < DataProcessingProcessML
    % Process Class for identifying quiescent windows by using Ljung-Box test.
    %
    % Updates:
    % J Noh, 7/2021. Modify QuiescentWindowDetectionProcessML for GCA pacakge. 
    % Qiongjing (Jenny) Zou, Aug 2018
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

    methods (Access = public)
        function obj = QuiescentWindowDetectionForGCAProcessML(owner, varargin)
            if nargin == 0
                super_args = {};
            else 
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);  % changed to ML's outputDir
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = QuiescentWindowDetectionForGCAProcessML.getName;
                super_args{3} = @QuiescentWindowDetectionForGCAMLWrapper;
                if isempty(funParams)
                    funParams=QuiescentWindowDetectionForGCAProcessML.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcessML(super_args{:});

        end        
    end

    methods (Static)
        function name = getName()
            name = 'Quiescent Window Detection';
        end

        function h = GUI()
            h= @QuiescentWindowDetectionForGCAProcessMLGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieList'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set up global parameters from SetUpGCAParametersProcessML
            setupProcID = owner.getProcessIndex('SetUpGCAParametersProcessML');
            funParams = parseProcessParams(owner.processes_{setupProcID});    

            funParams.movmeanNum = 5;   % # of windows for moving averages of P-values for smooth outcome
            funParams.topograph = 'off';
        end
    end
end