function varargout = QuiescentWindowDetectionForGCAProcessMLGUI(varargin)
%QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI MATLAB code file for QuiescentWindowDetectionForGCAProcessMLGUI.fig
%      QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI, by itself, creates a new QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI or raises the existing
%      singleton*.
%
%      H = QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI returns the handle to a new QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI or the handle to
%      the existing singleton*.
%
%      QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI('Property','Value',...) creates a new QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to QuiescentWindowDetectionForGCAProcessMLGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI('CALLBACK') and QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in QUIESCENTWINDOWDETECTIONFORGCAPROCESSMLGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Edit the above text to modify the response to help QuiescentWindowDetectionForGCAProcessMLGUI

% Last Modified by GUIDE v2.5 18-Jul-2021 18:56:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QuiescentWindowDetectionForGCAProcessMLGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @QuiescentWindowDetectionForGCAProcessMLGUI_OutputFcn, ...
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


% --- Executes just before QuiescentWindowDetectionForGCAProcessMLGUI is made visible.
function QuiescentWindowDetectionForGCAProcessMLGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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

% Disabled pars for User's information
set(handles.checkbox_lowFreqSubtraction, 'Enable', 'off')
set(handles.edit_avgLengthOfProtRetCycle, 'Enable', 'off')
set(handles.edit_chanCodeWithPreprocess1, 'Enable', 'off')
set(handles.edit_chanCodeWithPreprocess2, 'Enable', 'off')
set(handles.edit_chanNameWithPreprocess1, 'Enable', 'off')
set(handles.edit_chanNameWithPreprocess2, 'Enable', 'off')

if ~funParams.lowFreqSubtraction
  set(handles.checkbox_lowFreqSubtraction, 'Value', 0)
else
  set(handles.checkbox_lowFreqSubtraction, 'Value', 1)
end

set(handles.edit_avgLengthOfProtRetCycle, 'String', num2str(funParams.avgLengthOfProtRetCycle))
set(handles.edit_chanCodeWithPreprocess1, 'String', num2str(funParams.chanCodeWithPreprocess(1)))
set(handles.edit_chanCodeWithPreprocess2, 'String', num2str(funParams.chanCodeWithPreprocess(2)))
set(handles.edit_chanNameWithPreprocess1, 'String', funParams.chanNameWithPreprocess{1})
set(handles.edit_chanNameWithPreprocess2, 'String', funParams.chanNameWithPreprocess{2})

% 
if isequal(funParams.figFlag, 'on')
  set(handles.checkbox_figFlag, 'Value', 1)
else
  set(handles.checkbox_figFlag, 'Value', 0)
end

if isequal(funParams.topograph, 'on')
  set(handles.checkbox_topograph, 'Value', 1)
else
  set(handles.checkbox_topograph, 'Value', 0)
end

set(handles.edit_movmeanNum, 'String',num2str(funParams.movmeanNum)) 

% Choose default command line output for SNRoverSmoothingParamsProcessMLGUI
handles.output = hObject;

% Update handles structure
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);

% UIWAIT makes QuiescentWindowDetectionForGCAProcessMLGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = QuiescentWindowDetectionForGCAProcessMLGUI_OutputFcn(hObject, eventdata, handles)
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
if get(handles.checkbox_figFlag, 'Value') == 1
    funParams.figFlag = 'on';
else
    funParams.figFlag = 'off';
end

if get(handles.checkbox_topograph, 'Value') == 1
    funParams.topograph = 'on';
else
    funParams.topograph = 'off';
end

% edit text
funParams.movmeanNum = str2double(get(handles.edit_movmeanNum, 'String')); 

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);



function edit_movmeanNum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_movmeanNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_movmeanNum as text
%        str2double(get(hObject,'String')) returns contents of edit_movmeanNum as a double


% --- Executes during object creation, after setting all properties.
function edit_movmeanNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_movmeanNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function edit_chanCodeWithPreprocess1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chanCodeWithPreprocess1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chanCodeWithPreprocess1 as text
%        str2double(get(hObject,'String')) returns contents of edit_chanCodeWithPreprocess1 as a double


% --- Executes during object creation, after setting all properties.
function edit_chanCodeWithPreprocess1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chanCodeWithPreprocess1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_chanNameWithPreprocess1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chanNameWithPreprocess1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chanNameWithPreprocess1 as text
%        str2double(get(hObject,'String')) returns contents of edit_chanNameWithPreprocess1 as a double


% --- Executes during object creation, after setting all properties.
function edit_chanNameWithPreprocess1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chanNameWithPreprocess1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_chanNameWithPreprocess2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chanNameWithPreprocess2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chanNameWithPreprocess2 as text
%        str2double(get(hObject,'String')) returns contents of edit_chanNameWithPreprocess2 as a double


% --- Executes during object creation, after setting all properties.
function edit_chanNameWithPreprocess2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chanNameWithPreprocess2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_chanCodeWithPreprocess2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chanCodeWithPreprocess2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chanCodeWithPreprocess2 as text
%        str2double(get(hObject,'String')) returns contents of edit_chanCodeWithPreprocess2 as a double


% --- Executes during object creation, after setting all properties.
function edit_chanCodeWithPreprocess2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chanCodeWithPreprocess2 (see GCBO)
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


% --- Executes on button press in checkbox_lowFreqSubtraction.
function checkbox_lowFreqSubtraction_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_lowFreqSubtraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_lowFreqSubtraction


% --- Executes on button press in checkbox_topograph.
function checkbox_topograph_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_topograph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_topograph
