function varargout = SetUpGCAParametersProcessMLGUI(varargin)
%SETUPGCAPARAMETERSPROCESSMLGUI MATLAB code file for SetUpGCAParametersProcessMLGUI.fig
%      SETUPGCAPARAMETERSPROCESSMLGUI, by itself, creates a new SETUPGCAPARAMETERSPROCESSMLGUI or raises the existing
%      singleton*.
%
%      H = SETUPGCAPARAMETERSPROCESSMLGUI returns the handle to a new SETUPGCAPARAMETERSPROCESSMLGUI or the handle to
%      the existing singleton*.
%
%      SETUPGCAPARAMETERSPROCESSMLGUI('Property','Value',...) creates a new SETUPGCAPARAMETERSPROCESSMLGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SetUpGCAParametersProcessMLGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SETUPGCAPARAMETERSPROCESSMLGUI('CALLBACK') and SETUPGCAPARAMETERSPROCESSMLGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SETUPGCAPARAMETERSPROCESSMLGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Copyright (C) 2022, Danuser Lab - UTSouthwestern 
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

% Edit the above text to modify the response to help SetUpGCAParametersProcessMLGUI

% Last Modified by GUIDE v2.5 17-Jul-2021 15:07:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetUpGCAParametersProcessMLGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SetUpGCAParametersProcessMLGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SetUpGCAParametersProcessMLGUI is made visible.
function SetUpGCAParametersProcessMLGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

% Parameters setup 
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end          
funParams = userData.crtProc.funParams_;

if ~funParams.lowFreqSubtraction
  set(handles.checkbox_lowFreqSubtraction, 'Value', 0)
else
  set(handles.checkbox_lowFreqSubtraction, 'Value', 1)
end

if ~funParams.impute
  set(handles.checkbox_impute, 'Value', 0)
else
  set(handles.checkbox_impute, 'Value', 1)
end

if isequal(funParams.figFlag, 'on')
  set(handles.checkbox_figFlag, 'Value', 1)
else
  set(handles.checkbox_figFlag, 'Value', 0)
end

if ~funParams.WithN
  set(handles.checkbox_WithN, 'Value', 0)
else
  set(handles.checkbox_WithN, 'Value', 1)
end

set(handles.edit_maxLayer, 'String',num2str(funParams.maxLayer))
if isempty(funParams.avgLengthOfProtRetCycle) || isnan(funParams.avgLengthOfProtRetCycle)
    set(handles.edit_avgLengthOfProtRetCycle, 'String', [])
else
    set(handles.edit_avgLengthOfProtRetCycle, 'String',num2str(funParams.avgLengthOfProtRetCycle))
end

set(handles.edit_EWMAlambda, 'String',num2str(funParams.EWMAlambda))

% channel indexes, names
set(handles.edit_ChannelIndex1, 'String', num2str(funParams.ChannelIndex(1)))
set(handles.edit_ChannelIndex2, 'String', num2str(funParams.ChannelIndex(2)))
set(handles.edit_ChannelName1, 'String', funParams.ChannelName{1})
set(handles.edit_ChannelName2, 'String', funParams.ChannelName{2})

% Choose default command line output for SetUpGCAParametersProcessMLGUI
handles.output = hObject;

% Update handles structure
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);

% UIWAIT makes SetUpGCAParametersProcessMLGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SetUpGCAParametersProcessMLGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Delete figure
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

funParams = userData.crtProc.funParams_;

% -------- Process Sanity check --------
% ( only check underlying data )
try
    userData.crtProc.sanityCheck;
catch ME
    errordlg([ME.message 'Please double check your data'],...
                'Setting Error','modal');
    return;
end

% Retrieve GUI-defined parameters 
% checkboxes
funParams.lowFreqSubtraction = logical(get(handles.checkbox_lowFreqSubtraction, 'Value'));
funParams.impute = logical(get(handles.checkbox_impute, 'Value'));
if get(handles.checkbox_figFlag, 'Value') == 1
    funParams.figFlag = 'on';
else
    funParams.figFlag = 'off';
end
funParams.WithN = logical(get(handles.checkbox_WithN, 'Value'));

% edit text
funParams.maxLayer = str2double(get(handles.edit_maxLayer, 'String'));
funParams.avgLengthOfProtRetCycle = str2double(get(handles.edit_avgLengthOfProtRetCycle, 'String'));
funParams.EWMAlambda = str2double(get(handles.edit_EWMAlambda, 'String'));

% channel indexes, names
funParams.ChannelIndex(1) = str2double(get(handles.edit_ChannelIndex1, 'String'));
funParams.ChannelIndex(2) = str2double(get(handles.edit_ChannelIndex2, 'String'));
funParams.ChannelName{1} = get(handles.edit_ChannelName1, 'String');
funParams.ChannelName{2} = get(handles.edit_ChannelName2, 'String');

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% Run the SetUp process so that following processes can use the GCA global
% parameters 
userData.crtProc.run()



% --- Executes on button press in checkbox_impute.
function checkbox_impute_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_impute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_impute


% --- Executes on button press in checkbox_lowFreqSubtraction.
function checkbox_lowFreqSubtraction_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_lowFreqSubtraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_lowFreqSubtraction



function edit_maxLayer_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxLayer as text
%        str2double(get(hObject,'String')) returns contents of edit_maxLayer as a double


% --- Executes during object creation, after setting all properties.
function edit_maxLayer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_EWMAlambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_EWMAlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EWMAlambda as text
%        str2double(get(hObject,'String')) returns contents of edit_EWMAlambda as a double


% --- Executes during object creation, after setting all properties.
function edit_EWMAlambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_EWMAlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_figFlag.
function checkbox_figFlag_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_figFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_figFlag


% --- Executes on button press in checkbox_WithN.
function checkbox_WithN_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_WithN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_WithN



function edit_avgLengthOfProtRetCycle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_avgLengthOfProtRetCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_avgLengthOfProtRetCycle as text
%        str2double(get(hObject,'String')) returns contents of edit_avgLengthOfProtRetCycle as a double


% --- Executes during object creation, after setting all properties.
function edit_avgLengthOfProtRetCycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_avgLengthOfProtRetCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ChannelIndex1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ChannelIndex1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ChannelIndex1 as text
%        str2double(get(hObject,'String')) returns contents of edit_ChannelIndex1 as a double


% --- Executes during object creation, after setting all properties.
function edit_ChannelIndex1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ChannelIndex1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ChannelName1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ChannelName1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ChannelName1 as text
%        str2double(get(hObject,'String')) returns contents of edit_ChannelName1 as a double


% --- Executes during object creation, after setting all properties.
function edit_ChannelName1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ChannelName1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ChannelName2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ChannelName2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ChannelName2 as text
%        str2double(get(hObject,'String')) returns contents of edit_ChannelName2 as a double


% --- Executes during object creation, after setting all properties.
function edit_ChannelName2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ChannelName2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ChannelIndex2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ChannelIndex2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ChannelIndex2 as text
%        str2double(get(hObject,'String')) returns contents of edit_ChannelIndex2 as a double


% --- Executes during object creation, after setting all properties.
function edit_ChannelIndex2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ChannelIndex2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
