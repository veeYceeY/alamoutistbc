function varargout = alamouti_gui(varargin)
% ALAMOUTI_GUI MATLAB code for alamouti_gui.fig
%      ALAMOUTI_GUI, by itself, creates a new ALAMOUTI_GUI or raises the existing
%      singleton*.
%
%      H = ALAMOUTI_GUI returns the handle to a new ALAMOUTI_GUI or the handle to
%      the existing singleton*.
%
%      ALAMOUTI_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALAMOUTI_GUI.M with the given input arguments.
%
%      ALAMOUTI_GUI('Property','Value',...) creates a new ALAMOUTI_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alamouti_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alamouti_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alamouti_gui

% Last Modified by GUIDE v2.5 25-Apr-2015 19:17:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @alamouti_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @alamouti_gui_OutputFcn, ...
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


% --- Executes just before alamouti_gui is made visible.
function alamouti_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alamouti_gui (see VARARGIN)
global val_rx;
axes(handles.axes1);
img=imread('bg.jpg');
imshow(img);
val_rx=1;
% Choose default command line output for alamouti_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes alamouti_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = alamouti_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function N_box_Callback(hObject, eventdata, handles)
% hObject    handle to N_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_box as text
%        str2double(get(hObject,'String')) returns contents of N_box as a double
global N_bits;
N_bits = str2double(get(hObject,'string'));
if isnan(N_bits)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
	return
end

% --- Executes during object creation, after setting all properties.
function N_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global str_rx;
global val_rx;
global str_channel;
global val_channel;
global N_bits;
global gamma_alpha;
global gamma_gamma;
global log_mean;
global log_variance;
axes(handles.axes1);
img=imread('splash.png');
imshow(img);
pause(5);
nRx=val_rx;
%if (val_rx==1)
    switch str_channel{val_channel};
             case 'GammaGamma' 
                [simBer,Eb_N0_dB,theoryBer_nRx1,theoryBerMRC_nRx2,theoryBerAlamouti_nTx2_nRx1]=nRxscript_ber_2x2_alamouti_stbc_code_bpsk_gamma_gamma_channel(gamma_alpha,gamma_gamma,N_bits,nRx);
                axes(handles.axes1);
                title('BER for BPSK modulation with Alamouti STBC (gamma channel)');
             case 'LogNorm' 
                [simBer,Eb_N0_dB,theoryBer_nRx1,theoryBerMRC_nRx2,theoryBerAlamouti_nTx2_nRx1]=nRxscript_ber_2x2_alamouti_stbc_code_bpsk_logNorm_channel(log_mean,log_variance,N_bits,nRx);
                axes(handles.axes1);
                title('BER for BPSK modulation with Alamouti STBC (logNorm channel)');
             case 'Rayleigh' 
                [simBer,Eb_N0_dB,theoryBer_nRx1,theoryBerMRC_nRx2,theoryBerAlamouti_nTx2_nRx1]=nRxscript_ber_2x2_alamouti_stbc_code_bpsk_rayleigh_channel(N_bits,nRx);
                axes(handles.axes1);
                title('BER for BPSK modulation with Alamouti STBC (Rayleigh channel)');
    end



        semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
        hold on
        %semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2);
        semilogy(Eb_N0_dB,theoryBerAlamouti_nTx2_nRx1,'c+-','LineWidth',2);
        semilogy(Eb_N0_dB,simBer,'mo-','LineWidth',2);
        axis([0 25 10^-5 0.5])
        grid on
        legend('theory (nTx=1,nRx=1)','theory (nTx=2, nRx=1, Alamouti)', 'sim (nTx=2, nRx=1, Alamouti)');
        xlabel('Eb/No, dB');
        ylabel('Bit Error Rate');
        hold off

% --- Executes on selection change in channel_list.
function channel_list_Callback(hObject, eventdata, handles)
% hObject    handle to channel_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channel_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channel_list
global str_channel;
global val_channel;
str_channel = get(hObject, 'String');
val_channel = get(hObject,'Value');

% --- Executes during object creation, after setting all properties.
function channel_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mean_box_Callback(hObject, eventdata, handles)
% hObject    handle to Mean_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mean_box as text
%        str2double(get(hObject,'String')) returns contents of Mean_box as a double
global log_mean;
log_mean = str2double(get(hObject,'string'));
if isnan(log_mean)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
	return
end

% --- Executes during object creation, after setting all properties.
function Mean_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mean_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Variance_box_Callback(hObject, eventdata, handles)
% hObject    handle to Variance_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Variance_box as text
%        str2double(get(hObject,'String')) returns contents of Variance_box as a double
global log_variance;
log_variance = str2double(get(hObject,'string'));
if isnan(log_variance)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
	return
end

% --- Executes during object creation, after setting all properties.
function Variance_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Variance_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Alpha_box_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Alpha_box as text
%        str2double(get(hObject,'String')) returns contents of Alpha_box as a double
global gamma_alpha;
gamma_alpha = str2double(get(hObject,'string'));
if isnan(gamma_alpha)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
	return
end

% --- Executes during object creation, after setting all properties.
function Alpha_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gamma_box_Callback(hObject, eventdata, handles)
% hObject    handle to Gamma_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gamma_box as text
%        str2double(get(hObject,'String')) returns contents of Gamma_box as a double
global gamma_gamma;
gamma_gamma = str2double(get(hObject,'string'));
if isnan(gamma_gamma)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
	return
end

% --- Executes during object creation, after setting all properties.
function Gamma_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gamma_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in nRx_menu.
function nRx_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nRx_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nRx_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nRx_menu
global str_rx;
global val_rx;
str_rx = get(hObject, 'String');
val_rx = get(hObject,'Value');

% --- Executes during object creation, after setting all properties.
function nRx_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nRx_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
