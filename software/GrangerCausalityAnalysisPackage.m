classdef GrangerCausalityAnalysisPackage < Package
    % The main class of the GrangerCausalityAnalysis package
    %
    % Jungsik Noh, 7/2021.
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
        function obj = GrangerCausalityAnalysisPackage (owner,varargin)
            % Constructor of class GrangerCausalityAnalysisPackage
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar); % change to ML's outputDir
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                
                super_args{1} = owner;               
                super_args{2} = [outputDir filesep 'uIntegratePackage']; % Updated 2024-9-3. The old save folder name was GrangerCausalityAnalysisPackage. It was hard coded, need to change in all the places for this package.
            end
            % Call the superclass constructor
            obj = obj@Package(super_args{:}); 
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'u-integrate'; % Updated 2024-05-03. The old name was 'Granger-Causality Analysis'.
        end
        
        function m = getDependencyMatrix(i,j)
            %    1 2 3 4 5 6 7 8 9 {processes}
            m = [0 0 0 0 0 0 0 0 0;     % 1. SetUpGCAParametersProcessML
                 1 0 0 0 0 0 0 0 0;     % 2. SNRoverSmoothingParamsProcessML
                 1 0 0 0 0 0 0 0 0;     % 3. QuiescentWindowDetectionForGCAProcessML
                 1 0 2 0 0 0 0 0 0;     % 4. ActivityMapDescriptionForGCAProcessML
                 1 0 2 1 0 0 0 0 0;     % 5. XcorrAnalysisForGCAProcessML
                 1 0 2 0 0 0 0 0 0;     % 6. FluctuationProfilingForGCAProcessML
                 1 0 2 0 0 0 0 0 0;     % 7. GCA_3ComponentsProcessML
                 1 0 2 0 0 0 0 0 0;     % 8. GCA_2ComponentsProcessML
                 1 0 2 0 0 0 0 0 0];    % 9. GCA_1chInformationFlowProcessML
 
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = GrangerCausalityAnalysisPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            procContrs = {
                @SetUpGCAParametersProcessML, ...
                @SNRoverSmoothingParamsProcessML, ...
                @QuiescentWindowDetectionForGCAProcessML,...
                @ActivityMapDescriptionForGCAProcessML,...
                @XcorrAnalysisForGCAProcessML, ...
                @FluctuationProfilingForGCAProcessML, ...
                @GCA_3ComponentsProcessML, ...
                @GCA_2ComponentsProcessML, ...                
                @GCA_1chInformationFlowProcessML};
            
            if nargin==0, index=1:numel(procContrs); end
            procConstr=procContrs(index);
        end
        
        function classes = getProcessClassNames(index)
            classes = {
                'SetUpGCAParametersProcessML', ...
                'SNRoverSmoothingParamsProcessML', ...
                'QuiescentWindowDetectionForGCAProcessML',...
                'ActivityMapDescriptionForGCAProcessML',...
                'XcorrAnalysisForGCAProcessML',...
                'FluctuationProfilingForGCAProcessML', ...
                'GCA_3ComponentsProcessML', ...
                'GCA_2ComponentsProcessML', ...
                'GCA_1chInformationFlowProcessML'};
            if nargin==0, index=1:numel(classes); end
            classes=classes(index);
        end
        
        % add getMovieClass here, so will not call getMovieClass ('MovieData') in Package.m
        function class = getMovieClass() 
            % Retrieve the movie type on which the package can be applied
            class = 'MovieList';
        end
    end
end