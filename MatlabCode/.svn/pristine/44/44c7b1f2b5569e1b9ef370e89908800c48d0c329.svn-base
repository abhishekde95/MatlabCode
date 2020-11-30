%#ok<*DEFNU>

function IsoSampOnline()
global gl udpCom

KbName('UnifyKeyNames');

gl.modelout = []; % what the user wants to live in the NEX file as "model params"
gl.stimIdx = 1; % increments when a REQSTIMUPDATECD drops

header_recvd = false;

udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

% let the user pick the model and subject
[okay, gl.model_filepath] = user_select_model();
if ~okay, return; end
gl.model_revision = get_SVN_revision(gl.model_filepath);

% To generate the stimuli you must provide two things: 1) a MAT file with each subject's parameters
% and 2) a script that takes the subject's parameters and returns a matrix of stimuli (L M S TF) and
% a vector of model parameters you want sent to REX. Both the MAT file and the M file must have the
% same name (e.g., FOO.mat and FOO.m). It is up to you to ensure that this online program uses the
% correct files (i.e., keep MATLAB's path in mind). The `which` function is your friend.
% Specifically, you can use which('FOO', 'in', 'IsoSampOnline') to know exactly which 'FOO' is
% called from within this online program.

% The MAT file must have the following format: at the top level are structs named with the subject's
% identifer. The fields of the each subject's struct contain the relevant parameters to evaluate the
% model. The subject's struct is passed to the script in addition to arguments from REX. Here's an
% example layout of a MAT file and the corresponding script:

% >> s = load('FOO.mat')
% s =
%     Z: [1x1 struct]
%     A: [1x1 struct]
% >> s.Z
% ans =
%      model: [10x4 double]
%        sfs: [0.5 1 2 4]
%     domain: [3x2x4 double]
%
% Here is a brief description of the "modules": DTNT and LMTF
% function [stims,torex] = DTNT(target_nstims, subject_params, varargin)
%
% % `varargin` has additional arguments from REX - check the IsoSamp spot file
% % or private/LMTF.m
% <use the subject_params struct and REX arguments to compute the stimuli L M S TF>
% % you can follow {DTNT,LMTF}.m for examples of using isosurface + reducepatch
%
% stims = [lms tf]; torex = models(:,k);

gl.model = load_subjects_params(gl.model_filepath); % <-- 2-element cell array: 'LMTF' + structure with all models
if isempty(gl.model), return; end

codes = IsoSampCodes();

p = InitPStruct(0);

[udpCom.sock, success] = pnetStart(udpCom.port);
if ~success, return, end

server_handle = InitPlex();
message_REX('ONLINEINIT');

timetoleave = false;
while ~timetoleave
    timetoleave = user_quit() || deal_with_messages();
    [n, event_list] = PL_GetTS(server_handle);
    if n
        p = ProcessEventList(p, event_list);
        if ~header_recvd && ~isempty(find(p.events == codes.GAMMATABLECD, 1))
            parseHeader(codes, p);
            header_recvd = true;
        end
        p = parseTrialCodes(codes, p);
    end
    p = CleanUpEvents(p);
end
PL_Close(server_handle);

function yn = user_quit()
yn = false;
[keyisdown,~,keycodes] = KbCheck(-1);
if keyisdown
%    yn = all(keycodes([KbName('LeftControl') KbName('`~')]));
    yn = keycodes(KbName('Escape'));
end

function [okay,model_filepath] = user_select_model()
global gl

okay = true;
% reuse the previous model so users don't have to waste time clicking
if isfield(gl, 'model_filepath') && ~isempty(gl.model_filepath)
    model_filepath = gl.model_filepath;
else
    gl.model = {};
    
    [model_name, model_path] = ...
        uigetfile(fullfile(fileparts(mfilename('fullpath')), 'private', 'data', '*.mat'), ...
        'Which model parameters do you want to use in IsoSamp?');
    
    model_filepath = fullfile(model_path, model_name);
    if isequal(model_name,0) || isequal(model_path,0)
        okay = false;
        return
    end
end
[~,model_name] = fileparts(model_filepath);
fprintf('*** Loaded %s model parameters\n', model_name);

% Get the SVN revision number while maintaining that the module data and script are unmodified. This
% way can reproduce the state of the module offline.
function rev = get_SVN_revision(model_filepath)
[~,func_name] = fileparts(model_filepath);
script_filepath = which(func_name); % finds the first script with the same name as the MAT file
% execute an 'svn status' on the model's MAT and M files
% enforce that both files are in a "clean" state and return the working copy's revision number
% This function only returns the revision number associated with the *model params
% .mat file*, not the module (e.g. LMTF.m). The rationale here is that we
% expect the model params file to be updated periodically and we don't
% expect the module to change often if at all. GDLH and ZLB 4/12/16.
% (Previously, this returned the revision of the module script).
for file = {script_filepath model_filepath}
    status_cmd = sprintf('svn status -v "%s"', file{:});
    [rc,output] = system(status_cmd);
    if rc
        error('Had a problem executing ''%s''\n\toutput: %s\n', status_cmd, output);
    end
    output = regexp(output, '\n', 'split');
    output(cellfun('isempty', strfind(output, file{:}))) = []; % don't contain the filename? see ya
    if isempty(output), error('"%s" isn''t under version control', file{:}); end
    status_msg = output{1};
%    if status_msg(1) ~= ' ' % a blank first character means all is well
%        error('"%s" is in a non-clean state within version control; commit or revert!', file{:});
%    end
% GDLH Temporarily short-circuiting the above check since I'm actively
% working on LMTF.m (so it is "M"odified).
end
rev = str2double(regexp(status_msg, '\d+', 'match', 'once')); % the first consecutive set of digits

% Return the most recently used model parameters or show a dialog box to select the subject.
function model = load_subjects_params(model_filepath)
global gl

if isfield(gl, 'model') && ~isempty(gl.model)
    model = gl.model;
    gl.modelout = gl.model{2}.model;
    fprintf('*** Using %s model parameters for the most recently selected subject "%s".\n', model{1}, gl.subject_id);
    fprintf('*** Please close %s and execute ''clear all'' to select a different model or subject!\n\n', mfilename);
    return
end

s = load(model_filepath);
[~,script_name] = fileparts(model_filepath);
subject_ids = fieldnames(s);
[gl.subject,okay] = listdlg('Name', 'IsoSamp', 'PromptString', 'Which subject?', ...
    'SelectionMode', 'single', 'ListSize', [160 100], 'ListString', subject_ids);

if ~okay, model = []; return; end

gl.subject_id = subject_ids{gl.subject};
model_params = s.(gl.subject_id);
% GDLH added line below
% We can get model parameters early now
gl.modelout = model_params.model;
model = {script_name model_params};

function p = parseTrialCodes(codes, p)
global gl

% Was there a request to increment the stimulus index?
% This handshake is necessary to ensure stimIdx increments after a correct trial.
L = p.events == codes.REQSTIMUPDATECD;
if any(L)
    gl.stimIdx = gl.stimIdx + 1;
    codeidx = find(L, 1, 'last');
    p.lastprocessed_t = p.times(codeidx);
    message_REX('STIMIDXUPDATED');
end

function parseHeader(codes, p)
global gl

gl.bkgndrgb = GetValsIfPossible([codes.BKGNDRGBCD codes.BKGNDRGBCD], p, 'double');
spectra = reshape(GetValsIfPossible([codes.MONSPDCD codes.MONSPDCD], p, 'double'), [], 3);
fundamentals = reshape(GetValsIfPossible([codes.FUNDAMENTALSCD codes.FUNDAMENTALSCD], p, 'double'), [], 3);
P_device = SplineSpd(linspace(380,780,size(spectra,1))', spectra, ...
    linspace(380,780,size(fundamentals,1))');
gl.M = fundamentals'*P_device;
message_REX('GOTM_MATRIX');

% This function is the main entry point to Plexon from REX.
% The contents of varargin are from REX (check GENSTIMS_REQ and genStimSet() in IsoSamp.d)
function num_stims = genStimSet(target_nstims, use_prev_stims, varargin)
global gl

rf = varargin{5};
gl.stimIdx = 1;
model_str = gl.model{1};
subject_params = gl.model{2};
StimuliGenerator = str2func(model_str);
% The LMTF-StimuliGenerator function will reuse the model parameters if the RF matches!

% If the user wants to replay the previous set of stimuli, then make sure it's a sane choice.
if use_prev_stims
    if isfield(gl, 'prev_stims') && all(rf == gl.prev_rf)
        gl.stims = gl.prev_stims;
        gl.modelout = gl.prev_modelout;
        gl.localmodel = gl.prev_localmodel;
    elseif ~isfield(gl, 'prev_stims') % first run - no stimuli yet
        [gl.stims, gl.localmodel] = StimuliGenerator(target_nstims, subject_params, varargin{:});
    else
        warning('Can''t use previous IsoSamp stimuli because the RFs don''t match!\n');
        [gl.stims, gl.localmodel] = StimuliGenerator(target_nstims, subject_params, varargin{:});
    end
else
    [gl.stims gl.localmodel] = StimuliGenerator(target_nstims, subject_params, varargin{:});
end
num_stims = size(gl.stims, 1);
gl.stims = gl.stims(randperm(num_stims),:);
gl.prev_stims = gl.stims;
gl.prev_modelout = gl.modelout;
gl.prev_localmodel = gl.localmodel;
gl.prev_rf = rf;

% REX calls this; stimIdx updates when REX drops a REQSTIMUPDATECD. Check parseTrialCodes above.
function stim = nextStim()
global gl

if gl.stimIdx > size(gl.stims, 1)
    gl.stims = gl.stims(randperm(size(gl.stims, 1)),:);
    gl.stimIdx = 1;
end
stim = gl.stims(gl.stimIdx,:)

function leave = deal_with_messages()
global udpCom gl %#ok<NUSED>

leave = false;
bytes_avail = pnet(udpCom.sock, 'readpacket', 2000, 'noblock');
if bytes_avail
    message = pnet(udpCom.sock, 'read', bytes_avail, 'char');
    if strncmp(message, 'return', 6)
        stk = dbstack(); % Check whether called from another function or from command line
        if ~strcmp(stk(end).name, mfilename)
            leave = true;
        end
    end
    try
        eval(message);
    catch ex
        fprintf('Trouble with message: "%s"\n', message);
        disp(getReport(ex));
    end
end

function message_REX(message)
global udpCom
pnet(udpCom.sock, 'write', [message '>> >>']);
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
