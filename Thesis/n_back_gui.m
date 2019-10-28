function varargout = n_back_gui(varargin)
% N_BACK_GUI MATLAB code for n_back_gui.fig
%      N_BACK_GUI, by itself, creates a new N_BACK_GUI or raises the existing
%      singleton*.
%
%      H = N_BACK_GUI returns the handle to a new N_BACK_GUI or the handle to
%      the existing singleton*.
%
%      N_BACK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in N_BACK_GUI.M with the given input arguments.
%
%      N_BACK_GUI('Property','Value',...) creates a new N_BACK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before n_back_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to n_back_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help n_back_gui

% Last Modified by GUIDE v2.5 27-Oct-2019 17:01:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @n_back_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @n_back_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before n_back_gui is made visible.
function n_back_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to n_back_gui (see VARARGIN)


% initialize controller robot
import java.awt.Robot;
import java.awt.event.*;
handles.Markbot_mk_1.mouse = Robot;
handles.Markbot_mk_1.left_click = InputEvent.BUTTON1_MASK;
handles.Markbot_mk_1.right_click = InputEvent.BUTTON2_MASK;
handles.screenSize = get(0, 'screensize');

% launch STEM
%     STEM_controller(handles.Markbot_mk_1, "launch stem");

% Set default command line output for n_back_gui
global Start;
global Stop;
Start= false;
Stop= false;
handles.n_value= str2double(handles.n_value_popup.String{handles.n_value_popup.Value});
handles.n_correct= str2double(handles.n_correct_popup.String{handles.n_correct_popup.Value});
handles.deck_length= str2double(handles.length_popup.String{handles.length_popup.Value});
handles.output = hObject;

handles.mainscreen.FontSize= 20;
handles.mainscreen.HorizontalAlignment= 'center';
handles.mainscreen.String= "Welcome to the N-back working memory exercise! Select your test parameters from the bottom left corner of the screen, then press 'Start' to begin.";
handles.stimulus_number.String= "Stimulus number: ";
handles.your_answer.String= "Your answer: ";
guidata(hObject, handles);
% UIWAIT makes n_back_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = n_back_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in 'No Match".
function pushbutton_nomatch_Callback(hObject, eventdata, handles)
% hObject    handle to nomatch_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Display surf plot of the currently selected data.
global Start;

if Start== false
    current_text= handles.mainscreen.String;
    handles.mainscreen.String= "The exercise not started yet.";
    pause(1);
    handles.mainscreen.String= current_text;
else
    handles.your_answer.String= "Your answer: No Match";
end
guidata(hObject, handles);

% --- Executes on button press in "Match".
function pushbutton_match_Callback(hObject, eventdata, handles)
% hObject    handle to match_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Display mesh plot of the currently selected data.
global Start;

if Start== false
    current_text= handles.mainscreen.String;
    handles.mainscreen.String= "The exercise not started yet.";
    pause(1);
    handles.mainscreen.String= current_text;
else
    handles.your_answer.String= "Your answer: Match";
end
guidata(hObject, handles);

  
% --- Executes on button press in starts test.
function pushbutton_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Display contour plot of the currently selected data.
global Start;
global Stop;

if Start== false
%     Stop= false;
%     handles.mainscreen.FontSize= 20;
%     handles.mainscreen.HorizontalAlignment= 'center';
%     handles.mainscreen.String= "Welcome to the N-back working memory exercise! Select your test parameters from the bottom left corner of the screen, then press 'Start' to begin.";
%     handles.stimulus_number.String= "Stimulus number: ";
%     handles.your_answer.String= "Your answer: ";
%     Start= false;
% 
% % Start has been pressed once. Show the instructions
% elseif Start== false 
    [handles.n_back_deck, handles.answers]= build_nBack_deck(handles.n_value, handles.n_correct, handles.deck_length);
    handles.mainscreen.HorizontalAlignment= 'left';
    handles.mainscreen.String= ["Instructions:", newline,... 
                                "You will be shown a random series of letters. The goal of this exercise is to press 'Match' when the letter being shown is the same as the letter that was shown 'n' times ago. If the current letter does not match the letter shown 'n' times ago, press 'No Match' instead.", newline, "Press 'Start' again when are are ready to begin."];
    Start= true;

% Start has been pressed a second time. Start the exercise
else
    Stop= false;
    handles.mainscreen.FontSize= 20;
    handles.mainscreen.HorizontalAlignment= 'center';
    time= 5;
    while time> -1
        if Stop== true
            time= -1;
        elseif time== 1
            handles.mainscreen.String= "The exercise will begin in " + num2str(time) + " second.";
        else
            handles.mainscreen.String= "The exercise will begin in " + num2str(time) + " seconds.";
        end % end if
        pause(1);
        time= time-1;
    end % end while
    
    handles.mainscreen.String= " ";
    handles.mainscreen.FontSize= 400;
    
    for i= 1:handles.deck_length
        if Stop== false
            handles.stimulus_num= i;
            handles.stimulus_number.String= "Stimulus number: " + i + " of " + handles.deck_length;        
            STEM_controller(handles.Markbot_mk_1, "start capture");

            % present stimulus for 1 second
            handles.mainscreen.String= handles.n_back_deck(i);
            pause(1);

            % present fixation cross for 1 second
            handles.mainscreen.String= "+";
            pause(1);

    %         STEM_controller(handles.Markbot_mk_1, "stop capture", handles.answers(i), handles.your_answer.String);
        end % end if
    end % end for
    Stop= true;
    pushbutton_stop_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
end % end if



% --- Executes on button press in pauses test.
function pushbutton_pause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.pushbutton_pause.String== "Pause"
    handles.pushbutton_pause.String= "Unpause";
    pause();
else
    handles.pushbutton_pause.String= "Pause";
end % end if
guidata(hObject, handles);

% --- Executes on button press in ends test.
function pushbutton_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)  
global Start;
global Stop;

handles.mainscreen.FontSize= 20;
handles.mainscreen.HorizontalAlignment= 'center';

if Start== true % if it was started, stop it
    Stop= true;
    handles.mainscreen.String= "End of Exercise.";
    Start= false;
    pause(1.5);
elseif Start== false && Stop== true % if it wasn't started, don't stop it
    handles.mainscreen.String= "Welcome to the N-back working memory exercise! Select your test parameters from the bottom left corner of the screen, then press 'Start' to begin.";
    Stop= false;
elseif Start== false && Stop== false
    handles.mainscreen.String= "The exercise has not started yet.";
    Stop= false;
    pause(1.5);
% else
%     message= "Welcome to the N-back working memory exercise! Select your test parameters from the bottom left corner of the screen, then press 'Start' to begin.";
end % end if

handles.mainscreen.String= "Welcome to the N-back working memory exercise! Select your test parameters from the bottom left corner of the screen, then press 'Start' to begin.";
handles.stimulus_number.String= "Stimulus number: ";
handles.your_answer.String= "Your answer: ";
guidata(hObject, handles);


% --- Executes on selection change in selects value of n for test.
function popupmenu_nvalue_Callback(hObject, eventdata, handles)
% hObject    handle to n_value_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Determine the selected value.
str = get(hObject, 'String');

% Set the value of n
handles.n_value = str2double(str{get(hObject,'Value')});

% Save the handles structure.
guidata(hObject, handles)

% --- Executes on selection change in n back set length.
function length_popup_Callback(hObject, eventdata, handles)
% hObject    handle to length_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns length_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from length_popup
str = get(hObject, 'String');
handles.deck_length = str2double(str{get(hObject,'Value')});

if handles.deck_length <= handles.n_correct
    handles.n_correct= handles.deck_length- 5;
    handles.n_correct_popup.Value= find(strcmp(handles.n_correct_popup.String, num2str(handles.n_correct)));
end % end if
% Save the handles structure.
guidata(hObject, handles)

% --- Executes on selection change in n_correct_popup.
function n_correct_popup_Callback(hObject, eventdata, handles)
% hObject    handle to n_correct_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns n_correct_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n_correct_popup
str = get(hObject, 'String');
ncor = str2double(str{get(hObject,'Value')});

if handles.deck_length <= ncor
    handles.n_correct= handles.deck_length- 5;
    handles.n_correct_popup.Value= find(strcmp(handles.n_correct_popup.String, num2str(handles.n_correct)));
else
    handles.n_correct= ncor;
end % end if

% Save the handles structure.
guidata(hObject, handles)


% -------------------------------------------------------------------------
% Create Functions
% -------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function popupmenu_nvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_value_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function length_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function n_correct_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_correct_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton_match_CreateFcn(hObject, eventdata, handles)
% hObject    handle to match_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton_nomatch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nomatch_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton_pause_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes during object deletion, before destroying properties.
function axes1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
