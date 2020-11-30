% This function returns true if the machine is hooked up to a *Pixx device. The
% result is stored for fast retrieval on subsequent calls. Run `clear all` to
% remove the stored result.

function tf = IsVPixx()
persistent rslt
if isempty(rslt)
    if ismac
        cmd = 'system_profiler SPDisplaysDataType';
    elseif isunix % use an external script to get the display's name from the EDID
        this_path = fileparts(mfilename('fullpath'));
        this_path = fullfile(this_path, 'monitor-name.sh');
        cmd = ['bash ' this_path];
    else
        error('unknown operating system!');
    end
    
    [rc,response] = system(cmd);
    if rc ~= 0
        rslt = 0;
        warning(['Assuming this machine is NOT connected to a VPixx device!\n' ...
            'There was an error trying to figure out the connected display''s name.\n' ...
            'Here''s what was returned from the command:\n%s\n'], response);
    else % command successful -- is there a VPixx device attached?
        rslt = ~isempty(regexpi(response, 'pixx'));
    end
end
tf = rslt;
