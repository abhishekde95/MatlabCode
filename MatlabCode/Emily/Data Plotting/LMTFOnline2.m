function LMTFOnline2()
% If you're here wondering how LMTF works or to debug it, then the main meat of the
% online side is found in private/init_LMTF.m. There you can follow the comments
% and expand on any subfunction to go down the rabbit hole.
global gl udpCom

KbName('UnifyKeyNames');

header_recvd = false;
gl.init_new_exp = false;
gl.NSTIMPERROUND = 4;
gl.stimuli = [];
gl.HFIG = 111;

udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

codes = LMTFCodes();

p = InitPStruct(0);

udpCom.sock = pnetStart(udpCom.port);

if ~isfield(gl, 'data')
    gl.data = [];
    gl.isostim = [];
    gl.nfiles = 0;
    gl.processed_fnames = [];
    gl.isosamp_rf = [];
end

if ~isfield(gl, 'prev_path')
    gl.prev_path = nexfilepath('nexfilelists','Greg','LMTF');
end

[paradigmName, RFX, RFY, subjectID, notes] = DBInfoPopGUI;
fnames = fnamesFromDB(paradigmName, RFX, RFY, subjectID, notes);
if strcmp(fnames, 'cancel')
    gl.init_new_exp = true;
    gl.incoming_fnames = [];
    gl.std_inputs = false;
elseif strcmp(fnames, 'std_inputs')
    gl.init_new_exp = true;
    gl.incoming_fnames = [];
    gl.std_inputs = true;
else
    gl.incoming_fnames = findfile(fnames, 'N:\NexFiles');
    if ~isempty(gl.incoming_fnames)
        missing = cellfun('isempty', gl.incoming_fnames);
        if any(missing)
            missing_files = fnames(missing);
            warning('The following filenames were not found by findfile:');
            fprintf('\t%s\n', missing_files{:});
        end
    else
        warning('There were no NEX files in the loaded text file! Starting a new experiment...');
    end
end

server_handle = InitPlex();
c = onCleanup(@()PL_Close(server_handle));
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

function yn = user_quit()
yn = false;
[keyisdown,~,keycodes] = KbCheck(-1);
if keyisdown
    yn = keycodes(KbName('Escape'));
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
