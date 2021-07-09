function varargout = LBTestProcessMLGUI(varargin)
%LBTESTPROCESSMLGUI MATLAB code file for LBTestProcessMLGUI.fig
%      LBTESTPROCESSMLGUI, by itself, creates a new LBTESTPROCESSMLGUI or raises the existing
%      singleton*.
%
%      H = LBTESTPROCESSMLGUI returns the handle to a new LBTESTPROCESSMLGUI or the handle to
%      the existing singleton*.
%
%      LBTESTPROCESSMLGUI('Property','Value',...) creates a new LBTESTPROCESSMLGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to LBTestProcessMLGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LBTESTPROCESSMLGUI('CALLBACK') and LBTESTPROCESSMLGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LBTESTPROCESSMLGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LBTestProcessMLGUI

% Last Modified by GUIDE v2.5 28-Sep-2018 11:03:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LBTestProcessMLGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LBTestProcessMLGUI_OutputFcn, ...
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


% --- Executes just before LBTestProcessMLGUI is made visible.
function LBTestProcessMLGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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

if ~funParams.impute
  set(handles.checkbox_impute, 'Value', 0)
else
  set(handles.checkbox_impute, 'Value', 1)
end

if ~funParams.movingAvgSmoothing
  set(handles.checkbox_movingAvgSmoothing, 'Value', 0)
else
  set(handles.checkbox_movingAvgSmoothing, 'Value', 1)
end

if isequal(funParams.topograph, 'on')
  set(handles.checkbox_topograph, 'Value', 1)
else
  set(handles.checkbox_topograph, 'Value', 0)
end

if ~funParams.Folding
  set(handles.checkbox_folding, 'Value', 0)
else
  set(handles.checkbox_folding, 'Value', 1)
end

if ~funParams.derivative
  set(handles.checkbox_derivative, 'Value', 0)
else
  set(handles.checkbox_derivative, 'Value', 1)
end

set(handles.edit_smParam, 'String',num2str(funParams.smParam))

% Choose default command line output for LBTestProcessMLGUI
handles.output = hObject;

% Update handles structure
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = LBTestProcessMLGUI_OutputFcn(hObject, eventdata, handles)
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

% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end

% -------- Check user input --------
if isnan(str2double(get(handles.edit_smParam, 'String'))) ...
    || str2double(get(handles.edit_smParam, 'String')) < 0 ...
    || str2double(get(handles.edit_smParam, 'String')) > 1
  errordlg('Please provide a valid input for ''Activity map smoothing parameter for visualization''.','Setting Error','modal');
  return;
end

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
if get(handles.checkbox_impute, 'Value')
    funParams.impute = true;
else
    funParams.impute = false;
end

if get(handles.checkbox_movingAvgSmoothing, 'Value')
    funParams.movingAvgSmoothing = true;
else
    funParams.movingAvgSmoothing = false;
end

if get(handles.checkbox_topograph, 'Value')
    funParams.topograph = 'on';
else
    funParams.topograph = 'off';
end

if get(handles.checkbox_folding, 'Value')
    funParams.Folding = true;
else
    funParams.Folding = false;
end

if get(handles.checkbox_derivative, 'Value')
    funParams.derivative = true;
else
    funParams.derivative = false;
end

funParams.smParam = str2double(get(handles.edit_smParam, 'String'));

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


% --- Executes on button press in checkbox_movingAvgSmoothing.
function checkbox_movingAvgSmoothing_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_movingAvgSmoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_movingAvgSmoothing



function edit_smParam_Callback(hObject, eventdata, handles)
% hObject    handle to edit_smParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_smParam as text
%        str2double(get(hObject,'String')) returns contents of edit_smParam as a double


% --- Executes during object creation, after setting all properties.
function edit_smParam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_smParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_impute.
function checkbox_impute_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_impute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_impute


% --- Executes on button press in checkbox_topograph.
function checkbox_topograph_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_topograph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_topograph


% --- Executes on button press in checkbox_folding.
function checkbox_folding_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_folding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_folding


% --- Executes on button press in checkbox_derivative.
function checkbox_derivative_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_derivative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_derivative
