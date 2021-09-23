function varargout = XcorrAnalysisProcessMLGUI(varargin)
%XCORRANALYSISPROCESSMLGUI MATLAB code file for XcorrAnalysisProcessMLGUI.fig
%      XCORRANALYSISPROCESSMLGUI, by itself, creates a new XCORRANALYSISPROCESSMLGUI or raises the existing
%      singleton*.
%
%      H = XCORRANALYSISPROCESSMLGUI returns the handle to a new XCORRANALYSISPROCESSMLGUI or the handle to
%      the existing singleton*.
%
%      XCORRANALYSISPROCESSMLGUI('Property','Value',...) creates a new XCORRANALYSISPROCESSMLGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to XcorrAnalysisProcessMLGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      XCORRANALYSISPROCESSMLGUI('CALLBACK') and XCORRANALYSISPROCESSMLGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in XCORRANALYSISPROCESSMLGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help XcorrAnalysisProcessMLGUI

% Last Modified by GUIDE v2.5 06-Nov-2018 14:43:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @XcorrAnalysisProcessMLGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @XcorrAnalysisProcessMLGUI_OutputFcn, ...
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


% --- Executes just before XcorrAnalysisProcessMLGUI is made visible.
function XcorrAnalysisProcessMLGUI_OpeningFcn(hObject, eventdata, handles, varargin)
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

% Set up available input movieData if owner of the process is movieList (QZ added Nov 2018)
set(handles.listbox_availableMovieData,'String',userData.ML.movieDataFile_, ...
    'UserData',1:numel(userData.ML.movieDataFile_));

movieDataIndex = funParams.MovieDataIndex;

if ~funParams.impute
  set(handles.checkbox_impute, 'Value', 0)
else
  set(handles.checkbox_impute, 'Value', 1)
end

if ~funParams.omittedWindows
  set(handles.checkbox_omittedWindowsByLBTest, 'Value', 0)
else
  set(handles.checkbox_omittedWindowsByLBTest, 'Value', 1)
end

if isequal(funParams.topograph, 'on')
  set(handles.checkbox_topograph, 'Value', 1)
else
  set(handles.checkbox_topograph, 'Value', 0)
end

if ~funParams.WithN
  set(handles.checkbox_withN, 'Value', 0)
else
  set(handles.checkbox_withN, 'Value', 1)
end

if ~funParams.Folding
  set(handles.checkbox_folding, 'Value', 0)
else
  set(handles.checkbox_folding, 'Value', 1)
end

set(handles.edit_maxLayer, 'String',num2str(funParams.maxLayer))
set(handles.edit_numPerm, 'String',num2str(funParams.numPerm))
if funParams.lagMax == 5
    set(handles.edit_lagMax, 'String', 'Data driven default')
else
    set(handles.edit_lagMax, 'String',num2str(funParams.lagMax))
end

allChanName = join(funParams.chanName);
set(handles.edit_chanNames, 'String',allChanName{1})

% Choose default command line output for XcorrAnalysisProcessMLGUI
handles.output = hObject;

% Update handles structure
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = XcorrAnalysisProcessMLGUI_OutputFcn(hObject, eventdata, handles)
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
if isempty(get(handles.listbox_selectedMovieData, 'String'))
    errordlg('Please select at least one input MovieData from ''Available MovieData''.','Setting Error','modal')
    return;
end

if isnan(str2double(get(handles.edit_maxLayer, 'String'))) ...
    || str2double(get(handles.edit_maxLayer, 'String')) < 0
  errordlg('Please provide a valid input for ''Analyze layers from 1 to ''.','Setting Error','modal');
  return;
end 

if isnan(str2double(get(handles.edit_numPerm, 'String'))) ...
    || str2double(get(handles.edit_numPerm, 'String')) < 0
  errordlg('Please provide a valid input for ''Number of permutations''.','Setting Error','modal');
  return;
end 

if isnan(str2double(get(handles.edit_lagMax, 'String'))) ...
    && ~isequal(get(handles.edit_lagMax, 'String'), 'Data driven default') ...
    || str2double(get(handles.edit_lagMax, 'String')) < 0
  errordlg('Please provide a valid input for ''Maximum correlation lag''.','Setting Error','modal');
  return;
end 

if isempty(get(handles.edit_chanNames, 'String'))
  errordlg('Please provide a valid input for ''Channel Name(s)''.','Setting Error','modal');
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
movieDataIndex = get(handles.listbox_selectedMovieData, 'Userdata');
funParams.MovieDataIndex = movieDataIndex;

if get(handles.checkbox_impute, 'Value')
    funParams.impute = true;
else
    funParams.impute = false;
end

if get(handles.checkbox_omittedWindowsByLBTest, 'Value')
    funParams.omittedWindows = true;
else
    funParams.omittedWindows = false;
end

if get(handles.checkbox_topograph, 'Value')
    funParams.topograph = 'on';
else
    funParams.topograph = 'off';
end

if get(handles.checkbox_withN, 'Value')
    funParams.WithN = true;
else
    funParams.WithN = false;
end

if get(handles.checkbox_folding, 'Value')
    funParams.Folding = true;
else
    funParams.Folding = false;
end

funParams.maxLayer = str2double(get(handles.edit_maxLayer, 'String'));
funParams.numPerm = str2double(get(handles.edit_numPerm, 'String'));
if isequal(get(handles.edit_lagMax, 'String'), 'Data driven default')
  funParams.lagMax = 5;
else 
  funParams.lagMax = str2double(get(handles.edit_lagMax, 'String'));
end
funParams.chanName = strsplit(get(handles.edit_chanNames, 'String'));

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


function edit_numPerm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numPerm as text
%        str2double(get(hObject,'String')) returns contents of edit_numPerm as a double


% --- Executes during object creation, after setting all properties.
function edit_numPerm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numPerm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in checkbox_impute.
function checkbox_impute_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_impute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_impute


% --- Executes on button press in checkbox_omittedWindowsByLBTest.
function checkbox_omittedWindowsByLBTest_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_omittedWindowsByLBTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_omittedWindowsByLBTest



function edit_chanNames_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chanNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chanNames as text
%        str2double(get(hObject,'String')) returns contents of edit_chanNames as a double


% --- Executes during object creation, after setting all properties.
function edit_chanNames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chanNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_topograph.
function checkbox_topograph_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_topograph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_topograph


% --- Executes on button press in checkbox_withN.
function checkbox_withN_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_withN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_withN


% --- Executes on button press in checkbox_folding.
function checkbox_folding_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_folding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_folding



function edit_lagMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lagMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lagMax as text
%        str2double(get(hObject,'String')) returns contents of edit_lagMax as a double


% --- Executes during object creation, after setting all properties.
function edit_lagMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lagMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_applytoall.
function checkbox_applytoall_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_applytoall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_applytoall


% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
contents1 = get(handles.listbox_availableMovieData, 'String');

mdIndex1 = get(handles.listbox_availableMovieData, 'Userdata');
mdIndex2 = get(handles.listbox_selectedMovieData, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedMovieData, 'String', contents1);
        mdIndex2 = mdIndex1;
    case 0
        set(handles.listbox_selectedMovieData, 'String', {}, 'Value',1);
        mdIndex2 = [ ];
end
set(handles.listbox_selectedMovieData, 'UserData', mdIndex2);
% update_data(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button

contents1 = get(handles.listbox_availableMovieData, 'String');
contents2 = get(handles.listbox_selectedMovieData, 'String');
id = get(handles.listbox_availableMovieData, 'Value');

% If channel has already been added, return;
mdIndex1 = get(handles.listbox_availableMovieData, 'Userdata');
mdIndex2 = get(handles.listbox_selectedMovieData, 'Userdata');

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        mdIndex2 = cat(2, mdIndex2, mdIndex1(i));
    end
end

set(handles.listbox_selectedMovieData, 'String', contents2, 'Userdata', mdIndex2);
% update_data(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_selectedMovieData,'String');
id = get(handles.listbox_selectedMovieData,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
mdIndex2 = get(handles.listbox_selectedMovieData, 'Userdata');
mdIndex2(id) = [ ];
set(handles.listbox_selectedMovieData, 'Userdata', mdIndex2);

% Refresh listbox
set(handles.listbox_selectedMovieData,'String',contents);