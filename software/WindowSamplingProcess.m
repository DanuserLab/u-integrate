classdef WindowSamplingProcess < ImageSamplingProcess
    %Process
    %
    % Hunter Elliott
    % 7/2010
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
    
    methods (Access = public)
        
        function obj = WindowSamplingProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = WindowSamplingProcess.getName;
                super_args{3} = @sampleMovieWindows;
                if isempty(funParams)
                    funParams=WindowSamplingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                
            end
            
            obj = obj@ImageSamplingProcess(super_args{:});
        end
        

    end
    methods (Static)
        function name =getName()
            name = 'Window Sampling';
        end
        function name= GUI()
            name =@windowSamplingProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);%Default is to sample all channels
            funParams.ProcessIndex = [];%Default is to use raw images
            funParams.SegProcessIndex = [];%Default is to use masks which were used in windowing.
            funParams.MaskChannelIndex = [];%Default is to use channel which was used for windowing.
            funParams.OutputName = '';%Default is to use raw images
            funParams.OutputDirectory = [outputDir  filesep 'window_sampling'];
            funParams.BatchMode = false;
        end
        function samplableInput = getSamplableInput()
            % List process output that can be sampled
            processNames = horzcat('Raw images','DoubleProcessingProcess',...
                repmat({'KineticAnalysisProcess'},1,3),'FlowAnalysisProcess');
            samplableOutput = {'','','netMap','polyMap','depolyMap','speedMap'};
            samplableInput=cell2struct(vertcat(processNames,samplableOutput),...
                {'processName','samplableOutput'});  
        end
        
        
    end
end