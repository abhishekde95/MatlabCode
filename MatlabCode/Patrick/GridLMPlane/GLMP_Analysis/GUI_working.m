function varargout = GUI_working(varargin)
% GUI_WORKING MATLAB code for GUI_working.fig
%      GUI_WORKING, by itself, creates a new GUI_WORKING or raises the existing
%      singleton*.
%
%      H = GUI_WORKING returns the handle to a new GUI_WORKING or the handle to
%      the existing singleton*.
%
%      GUI_WORKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_WORKING.M with the given input arguments.
%
%      GUI_WORKING('Property','Value',...) creates a new GUI_WORKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_working_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_working_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_working

% Last Modified by GUIDE v2.5 09-Nov-2011 11:36:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_working_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_working_OutputFcn, ...
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

% --- Executes just before GUI_working is made visible.
function GUI_working_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_working (see VARARGIN)

% Choose default command line output for GUI_working
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using GUI_working.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes GUI_working wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_working_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in plot1_pushbutton1.
function plot1_pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to plot1_pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
    case 6
        [allspikes, bincenters, LpMspikes, LmMspikes] = gethistdata;
        subplot(3,1,1); hold on; grid on;% axis([0 max(bincenters) 0 1]); axis 'auto y';
        title('All Responses')
        hist(allspikes,bincenters)
        subplot(3,1,2); hold on; grid on;% axis([0 max(bincenters) 0 1]); axis 'auto y';
        title('Responses to L + M Stimuli')
        ylabel('Number of Spikes')
        hist(LpMspikes,bincenters)
        subplot(3,1,3); hold on; grid on;% axis([0 max(bincenters) 0 1]); axis 'auto y';
        title('Responses to L - M Stimuli')
        hist(LmMspikes,bincenters)
        xlabel('Miliseconds from Stimulus Onset')
        %axes(handles.axes1);
end



% ----------------------------------------------------
function [allspikes, bincenters, LpMspikes, LmMspikes] = gethistdata()
% Inserting data for analysis

rawdata = nex2stro('Users/jpatrickweller/MATLAB/GridLMPlane/Datafiles/S102811002.nex'); % Symetric Trench

% Preallocate Space
trial.nspikes = nan(size(rawdata.trial,1),1);
trial.tspikes = rawdata.ras(:,1);
trial.stimon = rawdata.trial(:,4);
trial.stimoff = rawdata.trial(:,3);
trial.duration = rawdata.trial(:,3) - rawdata.trial(:,4);
trial.stim = [rawdata.trial(:,5) rawdata.trial(:,6)];
for n = 1:size(rawdata.trial,1)
    trial.nspikes(n) = sum(logical(rawdata.ras{n,1} > rawdata.trial(n,4)) & logical(rawdata.ras{n,1} < rawdata.trial(n,3)));
end
trial.fr = trial.nspikes ./ trial.duration;

[coordinates,m,idx] = unique([rawdata.trial(:,5) rawdata.trial(:,6)],'rows');

for n = 1:size(coordinates,1)
    par(n).coordinates = coordinates(n,:);
    par(n).frs = trial.fr(idx == n);
    if numel(par(n).frs) < 3
        par(n).frs(numel(par(n).frs)+1:3,1) = nan;
    end
    par(n).meanfr = nanmean(par(n).frs);
    par(n).varfrs = nanvar(par(n).frs);
end


ptMat = cat(1,par(:).coordinates);
%frMat = cat(1,par(:).meanfr);

transMat = [.05 .5; -.05 .5]; %This matrix transforms from a normalized cone opponent space to an elongated cone opponent space with vertices [(.45,.55),(.55,.45),(-.45,-.55),(-.55,-.45)] calculated by GH and JPW
tPtMat = transMat\ptMat';
tPtMat = tPtMat';
trial.tstim = tPtMat;

for n = 1:size(rawdata.trial,1)
    trial.nspikes(n) = sum(logical(rawdata.ras{n,1} > rawdata.trial(n,4)) & logical(rawdata.ras{n,1} < rawdata.trial(n,3)));
end
trial.fr = trial.nspikes ./ trial.duration;

for n = 1:size(trial.stimoff,1)
    spiketimes = round((trial.tspikes{n} - trial.stimon(n)) * 1000);
    spiketimes = spiketimes(spiketimes > 0 & spiketimes < trial.duration(n)*1000);
    trial.normspiketimes{n} = spiketimes;
end

LpMidx = abs(trial.tstim(:,2)) >= abs(trial.tstim(:,1));
LmMidx = abs(trial.tstim(:,1)) >= abs(trial.tstim(:,2));

allspikes = cat(1,trial.normspiketimes{:});
LpMspikes = cat(1,trial.normspiketimes{LpMidx});
LmMspikes = cat(1,trial.normspiketimes{LmMidx});
bincenters = 0:5:max(allspikes);



% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)', 'histograms'});
