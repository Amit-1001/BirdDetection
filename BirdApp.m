function varargout = BirdApp(varargin)
%BIRDAPP MATLAB code file for BirdApp.fig
%      BIRDAPP, by itself, creates a new BIRDAPP or raises the existing
%      singleton*.
%
%      H = BIRDAPP returns the handle to a new BIRDAPP or the handle to
%      the existing singleton*.
%
%      BIRDAPP('Property','Value',...) creates a new BIRDAPP using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to BirdApp_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      BIRDAPP('CALLBACK') and BIRDAPP('CALLBACK',hObject,...) call the
%      local function named CALLBACK in BIRDAPP.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BirdApp

% Last Modified by GUIDE v2.5 05-Dec-2021 14:04:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BirdApp_OpeningFcn, ...
                   'gui_OutputFcn',  @BirdApp_OutputFcn, ...
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


% --- Executes just before BirdApp is made visible.
function BirdApp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for BirdApp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BirdApp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BirdApp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in train_extract.
function train_extract_Callback(hObject, eventdata, handles)
% hObject    handle to train_extract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp1 = get(handles.no_birds,'string');
load DATA.mat DIRPATH
temp2 = get(handles.train_samp,'string');
temp3 = get(handles.test_samp,'string');

numberB = str2double(temp1); % number of birds
numberTr = str2double(temp2); % number of training 
numberTe = str2double(temp3);



b =1
opn1 = get(handles.mfcc,'value');
opn2 = get(handles.timber,'value');
%save DATA.mat temp
%save 'DATA.mat' 

n = 1;
% For mfcc
if opn1 == 1
    for i=1:numberB
        for j =1:numberTr
            filename = sprintf('%d_%d.wav', i , j);
            filename=[DIRPATH '\' filename];

            [x, fs] = audioread(filename);%convert into samples x= total number of samples fs = sampling freq
            if (size(x,2) == 2) % converted into mono
                x= x(:,1);
            end
            a = miraudio(x); %MIRToolbox object

            mf=mfcc_common(x,fs,13); % first 13 component of mfcc
            mf = mf'; %transpose
            Training_feature(n,1:13) = mf(1,1:13); % storing feature of each mfcc component
            n = n+1;
        end
    end

    
    numberTr = numberTr+1;
    m =1;

    %this loop for testing feature
    for i=1:numberB
        for j =numberTr:numberTe
            filename = sprintf('%d_%d.wav', i , j);
            filename=[DIRPATH '\' filename];

            [x, fs] = audioread(filename);%convert into samples x= total number of samples fs = sampling freq
            if (size(x,2) == 2) % converted into mono
                x= x(:,1);
            end
            a = miraudio(x); %MIRToolbox object
            if opn1 == 1
                mf=mfcc_common(x,fs,13); % first 13 component of mfcc
                mf = mf'; %transpose
                Testing_feature(m,1:13) = mf(1,1:13); % storing feature of each mfcc component
                m = m+1;
            end
        end
    end
    x=10
      %Traingin feature and Testing feature of MFCC
    save 'training_feature.mat' 'Training_feature' 
    %save 'DATA.mat' 'Training_feature' 
    save 'testing_feature.mat' 'Testing_feature'
    b = 10;

    %putting training feature into excel file
    Features = 'Feature.xlsx';
    filenameXl = [DIRPATH '\' Features];
    b=10;

    xlswrite(filenameXl,Training_feature,'Traning_feature');
    b = 10;  
    
  
end



%for timber
x = 10;
if opn2 == 1
    n = 1;
   
    for i=1:numberB   
        for j =1:numberTr
            XY = 1;
            filename = sprintf('%d_%d.wav', i , j);
            filename=[DIRPATH '\' filename];
        
            [x, fs] = audioread(filename);%convert into samples x= total number of samples fs = sampling freq
            if (size(x,2) == 2) % converted into mono
                x= x(:,1);
            end
            a = miraudio(x); %MIRToolbox object
            
            ZCR = mirzerocross(a);  %ZCR keeps track of change of sign of signal
            Training_feature(n,XY) = mirgetdata(ZCR);
            XY = XY + 1;
        
            ROLL = mirrolloff(a);   %mirrolloff calculates the spectral roll-off in Hz.
            Training_feature(n,XY) = mirgetdata(ROLL);
            XY = XY + 1;
            
            BRIGHT = mirbrightness(a);  %calculates the spectral brightness, i.e. the amount
            Training_feature(n,XY) = mirgetdata(BRIGHT);
            XY = XY + 1;


            spectrum = mirspectrum(a);
            peaks = mirpeaks(spectrum);
            REGULARITY = mirregularity(peaks);     %calculates the irregularity of a spectrum, i.e.,
           %       the degree of variation of the successive peaks of the spectrum.
            Training_feature(n,XY) = mirgetdata(REGULARITY);
            XY = XY + 1;
            
            spectrum = mirspectrum(a);
            peaks = mirpeaks(spectrum);
            z = mirroughness(peaks);    %calculates the roughness, due to beating phenomenon between close frequency peaks. 
            ROUGH = mirgetdata(z);
            ROUGH = ROUGH';
            Training_feature(n,XY) = ROUGH; 
            XY = XY + 1;
            n = n+1;
        end
    end
    
    
    %Testing feature
    
    numberTr = numberTr+1;
     p = 1;
  
    for i=1:numberB   
        for j =numberTr:numberTe
            XY = 1;
            filename = sprintf('%d_%d.wav', i , j);
            filename=[DIRPATH '\' filename];
        
            [x, fs] = audioread(filename);%convert into samples x= total number of samples fs = sampling freq
            if (size(x,2) == 2) % converted into mono
                x= x(:,1);
            end
            a = miraudio(x); %MIRToolbox object
            
            ZCR = mirzerocross(a);  %ZCR keeps track of change of sign of signal
            Testing_feature(p,XY) = mirgetdata(ZCR);
            XY = XY + 1;
        
            ROLL = mirrolloff(a);   %mirrolloff calculates the spectral roll-off in Hz.
            Testing_feature(p,XY) = mirgetdata(ROLL);
            XY = XY + 1;
            
            BRIGHT = mirbrightness(a);  %calculates the spectral brightness, i.e. the amount
            Testing_feature(p,XY) = mirgetdata(BRIGHT);
            XY = XY + 1;


            spectrum = mirspectrum(a);
            peaks = mirpeaks(spectrum);
            REGULARITY = mirregularity(peaks);     %calculates the irregularity of a spectrum, i.e.,
           %       the degree of variation of the successive peaks of the spectrum.
            Testing_feature(p,XY) = mirgetdata(REGULARITY);
            XY = XY + 1;
            
            spectrum = mirspectrum(a);
            peaks = mirpeaks(spectrum);
            z = mirroughness(peaks);    %calculates the roughness, due to beating phenomenon between close frequency peaks. 
            ROUGH = mirgetdata(z);
            ROUGH = ROUGH';
           Testing_feature(p,XY) = ROUGH; 
            XY = XY + 1;
           p = p+1;
        end
    end 
    save 'training_feature.mat' 'Training_feature' 
    %save 'DATA.mat' 'Training_feature' 
    save 'testing_feature.mat' 'Testing_feature' 
end


var =10;




% --- Executes on button press in test_extract.
function test_extract_Callback(hObject, eventdata, handles)
% hObject    handle to test_extract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function training_features_Callback(hObject, eventdata, handles)
% hObject    handle to training_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of training_features as text
%        str2double(get(hObject,'String')) returns contents of training_features as a double


% --- Executes during object creation, after setting all properties.
function training_features_CreateFcn(hObject, eventdata, handles)
% hObject    handle to training_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function testing_features_Callback(hObject, eventdata, handles)
% hObject    handle to testing_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of testing_features as text
%        str2double(get(hObject,'String')) returns contents of testing_features as a double


% --- Executes during object creation, after setting all properties.
function testing_features_CreateFcn(hObject, eventdata, handles)
% hObject    handle to testing_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mfcc.
function mfcc_Callback(hObject, eventdata, handles)
% hObject    handle to mfcc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mfcc


% --- Executes on button press in timber.
function timber_Callback(hObject, eventdata, handles)
% hObject    handle to timber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timber


% --- Executes on button press in select_DB.
function select_DB_Callback(hObject, eventdata, handles)
% hObject    handle to select_DB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DIRPATH = uigetdir();
set(handles.path,'string',DIRPATH);
save 'DATA.mat' 'DIRPATH'


function no_birds_Callback(hObject, eventdata, handles)
% hObject    handle to no_birds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of no_birds as text
%        str2double(get(hObject,'String')) returns contents of no_birds as a double


% --- Executes during object creation, after setting all properties.
function no_birds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_birds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function train_samp_Callback(hObject, eventdata, handles)
% hObject    handle to train_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of train_samp as text
%        str2double(get(hObject,'String')) returns contents of train_samp as a double


% --- Executes during object creation, after setting all properties.
function train_samp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to train_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function test_samp_Callback(hObject, eventdata, handles)
% hObject    handle to test_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of test_samp as text
%        str2double(get(hObject,'String')) returns contents of test_samp as a double


% --- Executes during object creation, after setting all properties.
function test_samp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to test_samp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function path_Callback(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path as text
%        str2double(get(hObject,'String')) returns contents of path as a double


% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in generate_btn.
function generate_btn_Callback(hObject, eventdata, handles)
% hObject    handle to generate_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%load 'DATA.mat'

temp1 = get(handles.no_birds,'string');

temp2 = get(handles.train_samp,'string');
temp3 = get(handles.test_samp,'string');

numberB = str2double(temp1); % number of birds
numberTr = str2double(temp2); % number of training bird
numberTe = str2double(temp3); % number of testing bird



x = 10;
for i=0:numberB-1 % target generation
    for j=numberTr*i:numberTr*i+numberTr-1
        Target_Tr(i+1,j+1)=1;
    end
end
% % % %%
%%TARGET FOR TESTING
%%
Target_Ts=zeros(numberB,numberTe-numberTr);
numberTr=numberTe-numberTr;
for i=0:numberB-1
    for j=numberTr*i:numberTr*i+numberTr-1
        Target_Ts(i+1,j+1)=1;
    end
end

save 'DATA.mat' 'Target_Tr' 'Target_Ts'  
x =90




% --- Executes on button press in MSVM.
function MSVM_Callback(hObject, eventdata, handles)
% hObject    handle to MSVM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MSVM


% --- Executes on button press in FFBPNN.
function FFBPNN_Callback(hObject, eventdata, handles)
% hObject    handle to FFBPNN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FFBPNN


% --- Executes on button press in training_btn.
function training_btn_Callback(hObject, eventdata, handles)
% hObject    handle to training_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%empty neural network

temp1 = get(handles.no_birds,'string');

temp2 = get(handles.train_samp,'string');
temp3 = get(handles.test_samp,'string');

numberB = str2double(temp1); % number of birds
numberTr = str2double(temp2); % number of training bird
numberTe = str2double(temp3); %number of Testing bird

setdemorandstream(491218382)
op6 = get(handles.MSVM,'value');
op7 = get(handles.FFBPNN,'value'); 

      
    x = 100;
        if(op7==1)
              
                net = feedforwardnet([20,20,4]); % for 13 its 55% 
                
                % Now the network is ready to be trained. The samples are automatically
                % divided into training, validation and test sets. The training set is
                % used to teach the network. Training continues as long as the network
                % continues improving on the validation set. The test set provides a
                % completely independent measure of network accuracy.
                load 'training_feature.mat' 'Training_feature'
                load 'DATA.mat'  'Target_Tr'

                Training_feature
                Target_Tr
                %figure(10)
                %view(net)

                [net,tr] = train(net,Training_feature',Target_Tr); %taking transpose of Traning_feature
                
                
                
                
                save 'net.mat' 'net'
                save 'Data.mat' 'tr'
               
        end
    x = 100;
%hiden layer neuron = (2/3 * input+Out)
%%






% --- Executes on button press in Test_btn.
function Test_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Test_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

op6 = get(handles.MSVM,'value');
op7 = get(handles.FFBPNN,'value'); 

    if(op6==1)
    
    end
    
    if(op7==1)
           load 'testing_feature.mat' 'Testing_feature'
           load 'Data.mat' 'tr' 
           load 'net.mat' 'net'
           
           res = sim(net,Testing_feature');
           
           
            disp(res)
                numRows = size(res,1);
          
                res(res~=repmat(max(res),numRows,1)) = 0 % for finding max in each row

        
    end
  x = 101;
