classdef SetUpGCAParametersProcessML < DataProcessingProcessML
    % Process Class to set up global parameters of GCA pipeline.
    %
    %
    % Jungsik Noh, 7/2021
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
    
    %properties
    %    Property1
    %end
    
    methods (Access = public)
        
        function obj = SetUpGCAParametersProcessML(owner,varargin)
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
                super_args{2} = SetUpGCAParametersProcessML.getName;
                super_args{3} = @SetUpGCAParameters;
                if isempty(funParams)
                    funParams=SetUpGCAParametersProcessML.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end            
            obj = obj@DataProcessingProcessML(super_args{:});
        end
    end
    
    methods (Static)      
        
        function name = getName()
            name = 'Set Up GCA Parameters';
        end
        
        function h = GUI()
            h= @SetUpGCAParametersProcessMLGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieList'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set up global parameters
            numMovieData = numel(owner.movieDataFile_);             
            funParams.MovieDataIndex = 1:numMovieData;              % (1)
            funParams.ChannelIndex = [1, 2];
            funParams.ChannelName = {'channelName1', 'channelName2'};
            
            funParams.lowFreqSubtraction = true;             
            funParams.avgLengthOfProtRetCycle = [];                 
            funParams.maxLayer = 2; % for ichan=0, needs to be 1.   % (6)
            funParams.EWMAlambda = 0.5; % unnecessary here
                        
            % MD level params:
            % Used first MD 
            load(owner.movieDataFile_{1}) % load sample MD.            
            funParams.timeInterval = MD.timeInterval_;     

            % Did not register outputdir on MD level to this process:
            % funParams.OutputDirectory{1} = [outputDir  filesep 'Xcorr_Analysis' filesep mapDescriptivesDirName];

            % ML level params:
            % do not separate funParams.OutputDirectory here, otherwise folder 
            % icon on Control Panel of this package won't work.
            % complete output folders see ML.processes_{x}.outFilePaths_
            funParams.OutputDirectory = outputDir;                 
                                                
            funParams.impute = true;                                
            %funParams.movingAvgSmoothing = false;    % in GCA, do not use movingAvgSmoothing (3*3 patch).
            funParams.figFlag = 'off';                              % (11)
            funParams.WithN = false;         % if true, it uses an alternative windowSampling result
%                   which is obtained by sampleMovieWindowsWithN.m and includes number
%                   of pixels for each windows. Used for better sampling for Biosensor activity. 
%                   Default is false.
            
            % derived parameters
            funParams.chanCodeWithPreprocess = [];
            funParams.chanNameWithPreprocess = [];                  % (14)
            
            % For other various preprocessing methods in
            % mapOutlierImputation.m, directly specify the
            % chanCodeWithPreprocess1 and chanNameWithPreprocess1.
            % e.g., low-Frequency of signals, chanCodeWithPreprocess{1} = 3x
        end
    end
end

