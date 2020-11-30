%NEXFILEPATH   Return the mounted location of our NEX files along with any user-specified paths.
%Example usage: NEXFILEPATH('nexfilelists','Zack','File.txt')
%               NEXFILEPATH('Zack','Sedna')
% or simply     NEXFILEPATH()

% ZALB 2013/10/07

function filepath = nexfilepath(varargin)
filepath = '';

if nargin > 0
    initpath = nexfilepath();
    if ~isempty(initpath)
        filepath = fullfile(initpath, varargin{:});
    end
else
    server_ip = '128.95.153.12';
    win_cmd_str = ['net view ' server_ip];
    nix_cmd_str = ['mount | grep ' server_ip ' | sed -n ''s/.*on //;s/ %s.*//p'''];

    %#ok<*ASGLU>
    if ispc % assumes the PC mounted "NO BACKUP" as a drive
        [nil,partial_path] = system(win_cmd_str);
        colon_pos = find(partial_path == ':');
        partial_path = partial_path(colon_pos-1:colon_pos);
    elseif ismac % macs have a different `mount` format
        [nil,partial_path] = system(sprintf(nix_cmd_str, '('));
    elseif isunix
        [nil,partial_path] = system(sprintf(nix_cmd_str, 'type'));
    end

    if ~isempty(partial_path)
        filepath = [strtrim(partial_path) filesep 'NexFiles'];
    else
        warning('nexfilepath:notfound', ['Couldn''t find the mounted location of NO BACKUP.\n' ...
            'Make sure the server is powered on, that you are connected to the network, and that ' ...
            'the static IP address in %s.m is correct.'], mfilename);
    end
end
