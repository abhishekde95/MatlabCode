%necessary! populated by lmtf gui initial - used by LMTFOnline2
function filenames = fnamesFromDB(paradigmName, RFX, RFY, subjectID, notes)
if strcmp(paradigmName, 'cancel')
    filenames = 'cancel';
    return;
elseif strcmp(paradigmName, 'std_inputs')
    filenames = 'std_inputs';
    return;
end
conn = database('Nex_Paradigm_Sort','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
paramKeys = {'rfX', 'rfY', 'subjID', 'notes'};
paramVals = {RFX, RFY, subjectID, notes};
paramMap = containers.Map(paramKeys, paramVals);
klist = keys(paramMap);
unusableValues = {'*', 'null', 'cancel'};
for i = 1:length(paramMap)
    currentValue = values(paramMap, klist(i));
    currentValue = currentValue{:};
    if iscell(currentValue)
        currentValue = currentValue{:};
    elseif ~ischar(currentValue)
        currentValue = int2str(currentValue);
    end
    if ismember(currentValue, unusableValues)
        remove(paramMap, klist{i});
    else
        continue;
    end
end
usableKeys = keys(paramMap);
if isempty(usableKeys)
    paramQueryAdjusted = '';
else
    usableValues = values(paramMap, usableKeys);
    paramQuery = paramOrganizer(usableKeys, usableValues);
    paramQueryAdjusted = ['WHERE ', paramQuery];
end
fNQuery = sprintf('SELECT fileID FROM %s %s', paradigmName, paramQueryAdjusted);
filenames = fetch(conn, fNQuery);
close(conn);
    function paramQuery = paramOrganizer(keys, values)
        usefulParamOrder = [];
        q = char(39);
        for j = 1:length(keys)
            if strcmpi(keys{j}, 'notes')
                temp_val = values{j}';
                split_val = sprintf('''%s'',', temp_val{:});
                val = split_val(1:end-1);
                usefulParamOrder = [usefulParamOrder, '%', keys{j}, '%', val];
                continue;
            end
            if iscell(values{j})
                val = strcat(q, values{j}{:}, q);
            elseif ischar(values{j})
                val = strcat(q, values{j}, q);
            else
                val = int2str(values{j});
            end
            usefulParamOrder = [usefulParamOrder, '%', keys{j}, '%', val];
        end
        usefulParamCell_plusExtra = strsplit(usefulParamOrder, '%');
        usefulParamCell = usefulParamCell_plusExtra(2:end);
        isNotes = sum(ismember(usefulParamCell, 'notes'));
        if isNotes
            notesIDX = find(isNotes);
            note_string = sprintf('%s IN(%s)', usefulParamCell{notesIDX}, usefulParamCell{notesIDX+1});
            usefulParamCell(notesIDX:notesIDX+1) = [];
        end
        paramQuery_subs = sprintf(' AND %s = %s', usefulParamCell{:});
        if ~isempty(note_string)
            paramQuery = sprintf('%s%s', note_string, paramQuery_subs);
        else
           paramQuery = paramQuery_subs(6:end); 
        end
    end
end