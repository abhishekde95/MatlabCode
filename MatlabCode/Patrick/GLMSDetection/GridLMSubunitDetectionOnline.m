function GridLMSubunitDetectionOnline()
global gl udpCom

disp('In GridLMSubunitDetection')


KbName('UnifyKeyNames');
header_recvd = false;
gl.init_new_exp = false;
% gl.NSTIMPERROUND = 4;
% gl.stimuli = [];
% gl.HFIG = 111;

udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

codes = GLMSCodes(); % check the name of this file

p = InitPStruct(0);

[udpCom.sock, success] = pnetStart(udpCom.port);
if ~success, return, end


%User-Defined Variables
gl.nPres =5;
datafile = 'N102413002.nex'; % +L-M Chromatic (#43)


% Unpack datafile and set up experiment
disp('Loading data...')
%library = '/Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMSubunit/Datafiles/';
library = 'C:\Documents and Settings\JPatrickWeller\My Documents\Dropbox\Patrick\GLMS Data\';
keyboard
rawdata = nex2stro([char(library) datafile]);
disp('Organizing data...')
OrganizeRawGridLMSubunitDNData(rawdata);
disp('Setting up stimuli...')
SetUpStimuli()



% if ~isfield(gl, 'data')
%     gl.data = [];
%     gl.nfiles = 0;
%     gl.processed_fnames = [];
% end

% [filelist,filelistpath] = uigetfile(nexfilepath('nexfilelists','Greg','LMTF','*.txt'), ...
%     'Please select the file list or press Cancel to start a new experiment');
% if isequal(filelist,0)
%     gl.init_new_exp = true;
%     gl.incoming_fnames = [];
% else
%     fnames = flatten(fnamesFromTxt2(fullfile(filelistpath, filelist)));
%     gl.incoming_fnames = findfile(fnames);
%     missing = cellfun('isempty', gl.incoming_fnames);
%     if any(missing)
%         missing_files = fnames(missing);
%         warning('The following filenames were not found by findfile:');
%         fprintf('\t%s\n', missing_files{:});
%     end
% end
% 
% if isempty(gl.incoming_fnames)
%     gl.init_new_exp = true;
% end

server_handle = InitPlex();
message_REX('ONLINEINIT');

timetoleave = false;
while ~timetoleave
    timetoleave = user_quit() || deal_with_messages();
    [n, event_list] = PL_GetTS(server_handle);
    if n
        p = ProcessEventList(p, event_list);
        if ~header_recvd && ~isempty(find(p.events == codes.HEADERREADYCD, 1))
            parse_header(p, codes);
            header_recvd = true;
        end
        p = parse_trial_codes(p, codes);
    end
    p = CleanUpEvents(p);
end
PL_Close(server_handle);
end

function yn = user_quit()
yn = false;
[keyisdown,~,keycodes] = KbCheck(-1);
if keyisdown
    yn = all(keycodes([KbName('LeftControl') KbName('`~')]));
end
end


function leave = deal_with_messages()
global udpCom gl %#ok<NUSED>

leave = false;
bytes_avail = pnet(udpCom.sock, 'readpacket', 1000, 'noblock');
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
end


function SetUpStimuli()
global udpCom gl GLMP


end



function getStimParams()
global udpCom gl GLMP


end






