classdef SNRoverSmoothingParamsProcessML < DataProcessingProcessML
    % Process Class to compute Signal-to-Noise-Ratios over EWMA smoothing parameters 
    % to check optimal levels of smoothing. 
    %
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
    
    %properties
    %    Property1
    %end
    
    methods (Access = public)
        
        function obj = SNRoverSmoothingParamsProcessML(owner,varargin)
            % Construct an instance of this class
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar); 
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = SNRoverSmoothingParamsProcessML.getName;
                super_args{3} = @SNRoverSmoothingParamsMLWrapper;
                if isempty(funParams)
                    funParams=SNRoverSmoothingParamsProcessML.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end            
            obj = obj@DataProcessingProcessML(super_args{:});
        end
    end
    
    methods (Static)      
        
        function name = getName()
            name = 'Signal-to-Noise-Ratio Analysis';
        end
        
        function h = GUI()
            h= @SNRoverSmoothingParamsProcessMLGUI;
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
            
            % Order for auto-regressive models to compute SNR
            % default is 0.5*avgLengthOfProtRetCycle or 10 (frs) if
            % ProtRetCycle is empty
            if ~isempty(funParams.avgLengthOfProtRetCycle)
                funParams.ARorder = round(funParams.avgLengthOfProtRetCycle / 2);
            else
                funParams.ARorder = 10;     % frames
            end

            % set the frame size for moving median normalization
            if funParams.lowFreqSubtraction      
                funParams.movMedFrameSize = funParams.avgLengthOfProtRetCycle;
            else
                funParams.movMedFrameSize = [];
            end            
        end
    end
end

