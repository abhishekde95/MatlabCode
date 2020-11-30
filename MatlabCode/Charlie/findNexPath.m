function filepath = findNexPath(nexPaths, fname)
    %uses a pre-existing mapping of nexfile names to their absolute paths
    %to return the absolute path. Should speed things up when looking for
    %files over the network.
    
    [~, fname] = fileparts(fname); % in case the user passes an absolute path
    fname = [fname '.nex'];
    
    idx = find(strcmpi(nexPaths.names, fname), 2); % find at most 2 true entries
    
    %package the result.
    if ~isempty(idx)
        if length(idx) > 1
            filepath = fname;
            fprintf('File <%s> has none or more than one entry in NexFiles directory\n', fname);
        else
            filepath = strrep(nexPaths.paths{idx}, '$:$', filesep);
            filepath = fullfile(nexfilepath, filepath);
        end
    else
        fprintf('    ******  The experiment library needs to be updated!!!  ******\n');
        filepath = findfile(fname);
    end
end
