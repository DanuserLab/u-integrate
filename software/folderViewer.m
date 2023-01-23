function mainFig = folderViewer(ML, varargin)
%FOLDERVIEWER creates a graphical interface to open the folders of the
%analysis output of a MovieList.
% The analysis output are both at MovieData level (if any) and MovieList
% level.
% 
% h = folderViewer(ML, 'procId', 2)
% folderViewer(ML, 'procId', 2)
%
% This function is created to display (open results folders) of a
% MovieList's process (determined by the outFilePaths_ parameter).
% 
% Input 
%
%   ML - a MovieList
% 
%   procId - the specified index of a process for displaying the results.
%
%   Optional parameters in param/value pairs
%
%   folderIndex - Integer. For movie list input. If 0 display the movie list
%   and its analysis. If non-zero, set the index of the movie to be
%   displayed. Default: 0.
%
%   showProcTag - logical: displays process tags along with Process names
%                          for disambiguation.
%
% Output:
%   
%   mainFig - the handle of the main control interface
%
%
% Qiongjing (Jenny) Zou, Oct 2018
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


% Check input
ip = inputParser;
ip.addRequired('ML',@(x) isa(x,'MovieList'));
ip.addParameter('procId',[],@isnumeric);
ip.addParameter('folderIndex',0,@isscalar);
ip.addParameter('showProcTag',1,@islogical);
ip.parse(ML,varargin{:});

% Generate the main figure
mainFig=figure('Name','FolderViewer','Position',[0 0 200 200],...
    'NumberTitle','off','Tag','figure1','Toolbar','none','MenuBar','none',...
    'Color',get(0,'defaultUicontrolBackgroundColor'),'Resize','off',...
    'DeleteFcn', @(h,event) deleteViewer());
userData=get(mainFig,'UserData');
set(mainFig, 'UserData', userData);

% Read the MovieList and process index input 
userData.ML=ip.Results.ML;
userData.folderIndex=ip.Results.folderIndex;
userData.procId = ip.Results.procId;
if ~isempty(ip.Results.procId)
    procId = userData.ML.getProcessIndex(class(userData.ML.processes_{ip.Results.procId}));
else
    error("Must specify a ''procId'' and it must be positive integers.")
end

%
hPosition = 10;

%% Get image/overlay panel size and resize them
panelsLength = 1250; % qz changed from 500 to 1250 to make the whole FilderViewer window wider
panelsHeight = 0;

%% Create movie panel
moviePanel = uipanel(mainFig,...
    'Title','','BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
    'Units','pixels','Tag','uipanel_movie','BorderType','none');

% Create movie location edit box
uicontrol(moviePanel,'Style','text','Position',[10 hPosition 40 20],...
    'String','Movie','Tag','text_movie');

%%
% Create popupmenu if input is a MovieList, else  list the movie path

folderPaths = userData.ML.processes_{procId}.outFilePaths_;
folderIndex=0:numel(folderPaths);


uicontrol(moviePanel,'Style','popupmenu','Position',[60 hPosition panelsLength-110 20],...
    'String',vertcat('Please choose a folder to view results',folderPaths(:)),'UserData',folderIndex,...
    'Value',find(userData.folderIndex==folderIndex),...
    'HorizontalAlignment','left','BackgroundColor','white','Tag','popup_movie',...
    'Callback',@(h,event) switchFolder(h,guidata(h)));
if userData.folderIndex==0, set(findobj(moviePanel,'Tag','text_movie'),'String','List'); end
    

% Add help button
set(0,'CurrentFigure',mainFig)
hAxes = axes('Units','pixels','Position',[panelsLength-50 hPosition  20 20],...
    'Tag','axes_help', 'Parent', moviePanel);
icons = loadLCCBIcons();
Img = image(icons.questIconData);
set(hAxes, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'), 'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn, 'UserData', struct('class','movieViewer'));

% Add copyright
hPosition = hPosition+30;
uicontrol(moviePanel,'Style','text','Position',[10 hPosition panelsLength-100 20],...
    'String',getLCCBCopyright(),'Tag','text_copyright',...
    'HorizontalAlignment','left');

% Get overlay panel size
moviePanelSize = getPanelSize(moviePanel);
moviePanelHeight =moviePanelSize(2);


    %% Resize panels and figure
    sz=get(0,'ScreenSize');
    maxWidth = panelsLength+20;
    maxHeight = panelsHeight+moviePanelHeight;

set(mainFig,'Position',[sz(3)/50 sz(4)/3.8 maxWidth maxHeight]);
set(moviePanel,'Position',[10 panelsHeight+10 panelsLength moviePanelHeight]);
% Update handles structure and attach it to the main figure
handles = guihandles(mainFig);
guidata(handles.figure1, handles);

% Create redraw callbacks
% userData.getFigure = @(figName) getFigure(handles, figName);
set(handles.figure1,'UserData',userData);
end


function switchFolder(hObject,handles)
userData=get(handles.figure1,'UserData');
props=get(hObject,{'UserData','Value'});
if isequal(props{1}(props{2}), userData.folderIndex),return;end
if isempty(userData.procId)
    errordlg('Please choose a folder to view result from the dropdown options.')
else
    % Use the OS-specific command to open result in exploration window
    outputDir = userData.ML.processes_{userData.procId}.outFilePaths_{props{1}(props{2})};
    if ispc
        winopen(outputDir);
    elseif ismac
        system(sprintf('open %s',regexptranslate('escape',outputDir)));
    elseif isunix
        status = system(sprintf('xdg-open "%s"',regexptranslate('escape',outputDir)));
        % If a non-zero integer is returned, then display a message box
        if(status)
            msgbox(sprintf('Results can be found under %s',regexptranslate('escape',outputDir)));
        end
    else
        msgbox(sprintf('Results can be found under %s',regexptranslate('escape',outputDir)));
        % SB: Following command not working under Ubuntu (as well as gnome-open
        % & nautilus)
        % system(sprintf('xdg-open %s',regexptranslate('escape',outputDir)));
    end
end
end

function size = getPanelSize(hPanel)
if ~ishandle(hPanel), size=[0 0]; return; end
a=get(get(hPanel,'Children'),'Position');
P=vertcat(a{:});
size = [max(P(:,1)+P(:,3))+10 max(P(:,2)+P(:,4))+20];
end

function deleteViewer()

tags = {'viewerFig','optionsFig','graphFig'};
for i = 1:numel(tags)
    h = findobj(0,'-regexp','Tag', tags{i});
    if ~isempty(h), delete(h); end
end
end
