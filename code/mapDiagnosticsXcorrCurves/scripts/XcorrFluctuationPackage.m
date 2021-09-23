classdef XcorrFluctuationPackage < Package
    % The main class of the XcorrFluctuation package
    %
    % Qiongjing (Jenny) Zou, Sep 2018
    
    methods (Access = public)
        function obj = XcorrFluctuationPackage (owner,varargin)
            % Constructor of class XcorrFluctuationPackage
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
                super_args{2} = [outputDir filesep 'XcorrFluctuationPackage']; % but this outputDir is not used after run.
            end
            % Call the superclass constructor
            obj = obj@Package(super_args{:}); 
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'Fluctuation Analysis';
        end
        
        function m = getDependencyMatrix(i,j)
            %    1 2 3 4 {processes}
            m = [0 0 0 0; % QuiescentWindowDetectionProcessML
                 2 0 0 0; % ActivityMapDescriptionProcessML
                 2 1 0 0; % XcorrAnalysisProcessML
                 2 0 0 0]; % FluctuationProfilingProcessML
 
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = XcorrFluctuationPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            procContrs = {
                @QuiescentWindowDetectionProcessML,...
                @ActivityMapDescriptionProcessML,...
                @XcorrAnalysisProcessML, ...
                @FluctuationProfilingProcessML};
            
            if nargin==0, index=1:numel(procContrs); end
            procConstr=procContrs(index);
        end
        
        function classes = getProcessClassNames(index)
            classes = {
                'QuiescentWindowDetectionProcessML',...
                'ActivityMapDescriptionProcessML',...
                'XcorrAnalysisProcessML',...
                'FluctuationProfilingProcessML'};
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