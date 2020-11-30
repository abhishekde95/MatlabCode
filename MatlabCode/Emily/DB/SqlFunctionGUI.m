%Never ended up finishing this, mostly because no one really needed it.
%Nearly done, though, if anyone wants to come back to it.
function SqlFunctionGUI()
stored_queries = '';
SqlGUI(stored_queries)
end

function SqlGUI(stored_queries)
%recdate?
global action_choice DB_choice tables table_choice conn;
conditions_map = containers.Map('KeyType','char','ValueType','char');
secondary_conditions_map = containers.Map('KeyType','char','ValueType','char');
checkedQueries = ''; queries = stored_queries; conditions_string = ''; tot_columns = ''; primary_actions = ''; secondary_actions = '';
f = figure('Name', 'What SQL functions would you like to do?', 'NumberTitle', 'off');
set(f, 'Position', [600, 400, 700, 400]);
server_options = {'IBM Informix', 'Microsoft SQL Server 2005', 'MySQL', 'Oracle oci7 drivers', 'Oracle oci8 drivers', 'Oracle 10', 'Oracle Thin', 'PostgreSQL', 'Sybase'};
%GROUP 1: MAIN INFO AND ACTIONS
actions_options = {'SELECT', 'INSERT', 'UPDATE', 'DELETE'};
action_choice_text = uicontrol(f,'Style', 'Text', 'String', 'Function Type: ',  'Position', [20, 345, 80, 40], 'visible', 'off');
action_choice_menu = uicontrol(f, 'Style', 'popup', 'String', actions_options,  'Max', 1,'Callback', @getActionParams, 'Position', [100, 350, 90, 40], 'visible', 'off');
change_DB = uicontrol(f, 'String', 'Change DB', 'Callback', @DBselection, 'Position', [440, 10, 100, 40], 'visible', 'off');
% CONNECT TO A SERVER BEFORE STARTING
try ping(conn)
    if isopen(conn)
        onArray = [action_choice_text, action_choice_menu, change_DB];
        set(onArray, 'visible', 'on');
    else
        close(conn);
        connectToServer();
    end
catch
    connectToServer();
end
    function connectToServer(hcbo, eventStruct)
        Server_connect_text = uicontrol(f, 'Style', 'Text', 'String', 'Connect to a Server: ', 'Position', [20, 340, 150, 40]);
        Server_network_choice = uicontrol(f, 'Style', 'Edit', 'String', 'Enter the website or IP of your server', 'Position', [170, 350, 200, 40]);
        Server_Type_String = uicontrol(f, 'Style', 'Text', 'String', 'Select Server Type: ', 'Position', [20, 300, 150, 40]);
        Server_type = uicontrol(f, 'Style', 'popupmenu', 'String', server_options, 'Position', [170, 300, 90, 40]);
        Server_DB_string = uicontrol(f, 'Style', 'Text', 'String', 'Database Name: ', 'Position', [20, 260, 150, 40]);
        Server_DB_choice = uicontrol(f, 'Style', 'Edit', 'Position', [170, 270, 200, 40]);
        Server_UN_string = uicontrol(f, 'Style', 'Text', 'String', 'Username: ', 'Position', [20, 200, 150, 40]);
        Server_username = uicontrol(f, 'Style', 'Edit', 'String', '(leave blank if none)', 'Position', [170, 210, 200, 40]);
        Server_PW_string = uicontrol(f, 'Style', 'Text', 'String', 'Password: ', 'Position', [20, 150, 150, 40]);
        Server_password = uicontrol(f, 'Style', 'Edit', 'String', '(leave blank if none)', 'Position', [170, 160, 200, 40]);
        Connect_button = uicontrol(f, 'String', 'CONNECT', 'callback', 'uiresume(gcbf)', 'Position', [170, 50, 90, 40]);
        uiwait;
        network = get(Server_network_choice, 'String');
        type_IDX = get(Server_type, 'value');
        type = server_options{type_IDX};
        DB_choice = get(Server_DB_choice, 'String');
        UN = get(Server_username, 'String');
        PW = get(Server_password, 'String');
        if strcmp(UN, '(leave blank if none)')
            UN = '';
        end
        if strcmp(PW, '(leave blank if none)')
            PW = '';
        end
        conn = database(DB_choice, UN, PW, 'Vendor', type, 'Server', network);
        if ~isempty(conn.Message)
            if strfind(conn.Message, 'communications link failure')
                w1 = warndlg('Server Name Invalid', 'Try Again');
            elseif strfind(conn.Message, 'Vendor must be one of')
                w1 = warndlg('Incorrect Server Type', 'Try Again');
            else
                w1 = warndlg(conn.Message, 'Try Again');
            end
            waitfor(w1);
            AQ;
            return;
        end
        offArray = [Server_connect_text,Server_network_choice,Server_Type_String, Server_type,Server_DB_string,Server_DB_choice,Server_UN_string, Server_username,Server_PW_string,Server_password,Connect_button];
        set(offArray, 'visible', 'off');
        onArray = [action_choice_text, action_choice_menu, change_DB];
        set(onArray, 'visible', 'on');
    end

    function getActionParams(hcbo, eventStruct)
        actionIDX = get(action_choice_menu, 'value');
        action_choice = actions_options{actionIDX};
        onChoice();
    end
    function DBselection(hcbo, eventStruct)
        DBs = fetch(conn, 'SHOW DATABASES');
        DB_IDX = listdlg('ListString', DBs, 'SelectionMode', 'single','PromptString', 'Select another database from the server:');
        DB_choice = DBs{DB_IDX};
        AQ();
    end
    function onChoice(hcbo, eventStruct)
        tables = fetch(conn, sprintf('SHOW TABLES FROM %s', DB_choice));
        if isempty(tables)
            DBselection();
        end
        set(select_table, 'String', tables);
        offArray = [action_choice_menu, action_choice_text];
        set(offArray, 'visible', 'off');
        switch action_choice
            case 'SELECT' %SELECT <sel_opt> FROM <tables> <WHEREbtn>
                handlesArray = [Cancel, query_complete, Action_handle, select_table_cols, modifier_handle, select_table, WHERE_button];
                set(Action_handle, 'String', 'SELECT');
                set(select_table_cols, 'position', [70, 280, 120, 60]);
                set(modifier_handle, 'position', [190, 300, 50, 40]);
                set(select_table, 'position', [240, 280, 120, 60]);
                set(WHERE_button, 'position', [360, 280, 80, 40]);
            case 'INSERT' %INSERT INTO <tables> VALUES <condition> = <value>
                handlesArray = [Cancel, query_complete,Action_handle, select_table_cols, modifier_handle, select_table, equals_handle];
                set(Action_handle, 'String', 'INSERT INTO');
                set(Action_handle, 'position', [20, 300, 90, 40]);
                set(select_table, 'position', [100, 280, 120, 60]);
                set(modifier_handle, 'String', 'VALUES');
                set(modifier_handle, 'position', [220, 300, 50, 40]);
                set(select_table_cols, 'position', [270, 280, 120, 60]);
                set(equals_handle, 'position', [390, 280, 20, 40]);
            case 'UPDATE' %UPDATE <tables> SET <col_opt> <WHEREbtn>
                handlesArray = [Cancel, query_complete,Action_handle, select_table_cols, modifier_handle, select_table, equals_handle, WHERE_button];
                set(Action_handle, 'String', 'UPDATE');
                set(select_table, 'position', [70, 280, 120, 60]);
                set(modifier_handle, 'String', 'SET');
                set(modifier_handle, 'position', [190, 300, 30, 40]);
                set(select_table_cols, 'position', [220, 280, 120, 60]);
                set(equals_handle, 'position', [340, 280, 20, 40]);
                set(WHERE_button, 'position', [20, 200, 80, 40]);
            case 'DELETE' %DELETE FROM <tables> <WHEREbtn>
                handlesArray = [Cancel, query_complete,Action_handle, select_table, WHERE_button];
                set(Action_handle, 'String', 'DELETE FROM');
                set(Action_handle, 'position', [20, 300, 90, 40]);
                set(select_table, 'position', [110, 300, 120, 60]);
                set(WHERE_button, 'position', [230, 300, 80, 40]);
            otherwise
                warndlg('Invalid action choice');
        end
        set(handlesArray, 'visible', 'on');
    end
%INFO - GROUP 2 "ACTION OPTIONS"
Action_handle = uicontrol(f, 'Style', 'Text', 'String', 'SELECT', 'Position', [20, 300, 50, 40], 'visible', 'off');
modifier_handle = uicontrol(f, 'Style', 'Text', 'String', 'FROM', 'Position', [100, 300, 80, 40], 'visible', 'off');
select_table = uicontrol(f, 'Style', 'listbox', 'Max', 10, 'Callback', @colVal, 'Position', [280, 350, 90, 40], 'visible', 'off');
select_table_cols = uicontrol(f, 'Style', 'listbox', 'String', 'choose table first', 'Max', 100, 'Callback', @selectedChoice, 'Position', [180, 300, 90, 40], 'visible', 'off');
equals_handle = uicontrol(f, 'style', 'text', 'String', '=', 'visible', 'off');
another_equals = uicontrol(f, 'style', 'text', 'String', '=', 'visible', 'off');
WHERE_button = uicontrol(f, 'String', 'WHERE...', 'Position', [270, 200, 90, 40], 'Callback', @getMoreSpecs, 'visible', 'off');
Cancel = uicontrol(f, 'String', 'CANCEL', 'Callback', @AQ, 'Position', [150, 50, 50, 40], 'visible', 'off');
CloseAll = uicontrol(f, 'String', 'CLOSE AND DELETE ALL', 'Callback', 'close gcf', 'position', [540, 10, 150, 40]);
secondary_action = uicontrol(f, 'Style', 'listbox', 'String', tot_columns, 'Max', 100, 'Callback', @selectedChoice, 'visible', 'off');
query_complete = uicontrol(f, 'String', 'Query Complete', 'Position', [50, 50, 100, 40],'callback', @toQueryString, 'visible', 'off');
Drop = uicontrol(f, 'String', 'DROP SELECTED QUERIES', 'Callback', @drop, 'Position', [180, 10, 150, 40], 'visible', 'off');
Another = uicontrol(f, 'String', 'ADD ANOTHER QUERY', 'Callback', @AQ, 'Position', [20, 10, 150, 40], 'visible', 'off');
DB_send = uicontrol(f, 'String', 'SEND SELECTED QUERIES TO DB', 'callback', @done, 'Position', [330, 10, 200, 40], 'visible', 'off');
%GROUP 3: FUNCTIONS
    function colVal(hcbo, eventStruct)
        tableIDX = get(select_table, 'value');
        table_choice = (tables(tableIDX));
        table_choice_string = '';
        for q = 1:length(table_choice)
            table_choice_string = sprintf('%s,''%s''', table_choice_string, table_choice{q});
        end
        colQuery = sprintf('SELECT column_name FROM information_schema.columns WHERE table_schema = ''%s'' AND table_name IN(%s);', DB_choice, table_choice_string(2:end));
        tot_columns = fetch(conn, colQuery);
        tot_columns = vertcat('all', tot_columns);
        set(select_table_cols, 'string', tot_columns);
    end
    function selectedChoice(hcbo, eventStruct)
        if strcmpi(get(WHERE_button, 'String'), 'selected')
            optionIDX = get(secondary_action, 'value');
            if ~isempty(secondary_actions)
                secondary_actions = '';
            end
            secondary_actions = tot_columns(optionIDX);
            getValues(hcbo, eventStruct);
        else
            optionIDX = get(select_table_cols, 'value');
            if ~isempty(primary_actions)
                primary_actions = '';
            end
            primary_actions = tot_columns(optionIDX);
            if ~strcmp(action_choice, 'SELECT')
                getValues(hcbo, eventStruct);
            end
        end
    end
    function getValues(hcbo, eventStruct)
        if strcmpi(get(WHERE_button, 'String'), 'selected')
            final_pos = get(another_equals, 'position');
            items = secondary_actions;
        else
            final_pos = get(equals_handle, 'position');
            items = primary_actions;
        end
        if ismember(items, 'all')
            items = tot_columns(2:end);
        end
        x = final_pos(1)+30;
        y = final_pos(2)+40;
        existingBoxes = findall(f, 'Style', 'edit');
        delete(existingBoxes);
        for j = 1:length(items)
            editBox(j) = uicontrol(f,'style','edit','max',1,'string',items{j}, 'position', [x y 100 20]); %'keyPressFcn', {@boxContentCheck, j},
            y = y-20;
            if y<=50
                x=x+100;
                y= final_pos(2)+40;
            end
        end
        %         function boxContentCheck(hcbo, eventstruct, j)
        %             strtmp = get(editBox{j}, 'String');
        %             currType = items{j};
        %             switch currType
        %                 case 1
        %                     keyboard;
        %                 otherwise
        %                     keyboard;
        %             end
        %         end
        set(query_complete, 'Callback', @sortValues);
    end
    function sortValues(hcbo, eventStruct)
        editBoxes = findall(f, 'Style', 'edit');
        for k = 1:length(secondary_actions)
            conditions{k} = editBoxes(k).String;
        end
        conditions = fliplr(conditions);
        swapPoint = length(primary_actions);
        try
            if strcmp(action_choice, 'SELECT')
                dummy_conditions = repmat({''}, size(primary_actions));
                conditions_map = containers.Map(primary_actions, dummy_conditions);
            else
                conditions_map = containers.Map(primary_actions, conditions(1:swapPoint));
            end
        catch
            keyboard;
        end
        if ~isempty(secondary_actions) && strcmp(action_choice, 'SELECT')
            secondary_conditions_map = containers.Map(secondary_actions', conditions);
        elseif ~isempty(secondary_actions) && ~strcmp(action_choice, 'SELECT')
            secondary_conditions_map = containers.Map(secondary_actions', conditions(swapPoint+1:end));
        end
        toQueryString();
    end
    function toQueryString(hcbo, eventStruct)
        set(findall(f, 'type', 'uicontrol'), 'visible', 'off');
        onArray = [Drop, Another, DB_send, CloseAll];
        set(onArray, 'visible', 'on');
        if ~isempty(conditions_map)
            conditions_string = mapToString(conditions_map, 1);
        else
            conditions_string = '';
        end
        if ~isempty(secondary_conditions_map)
            secondary_conditions_string = mapToString(secondary_conditions_map, 1);
            if ~strcmp(action_choice, 'SELECT')
                conditions_string = [conditions_string(7:end), secondary_conditions_string];
            else
                conditions_string = secondary_conditions_string;
            end
        end
        selected_cols_string = cellToString(primary_actions);
        keyboard;
        tables_string = cellToString(table_choice);
        switch action_choice
            case 'SELECT'
                if strcmp(selected_cols_string, 'all')
                    queries{end+1} = sprintf('SELECT * FROM %s %s;', selected_cols_string, tables_string, conditions_string);
                else
                    queries{end+1} = sprintf('SELECT %s FROM %s %s;', selected_cols_string, tables_string, conditions_string);
                end
            case 'INSERT'
                queries{end+1} = sprintf('INSERT INTO %s VALUES %s;', tables_string, conditions_string(6:end));
            case 'UPDATE'
                queries{end+1} = sprintf('UPDATE %s SET %s;', tables_string, conditions_string);
            case 'DELETE'
                if strcmpi(get(WHERE_button, 'String'), 'selected')
                    queries{end+1} = sprintf('DELETE FROM %s %s;', tables_string, conditions_string);
                else
                    queries{end+1} = sprintf('DELETE %s;', tables_string);
                end
            otherwise
                warndlg('action choice error');
        end
        displayQueriesPanel();
    end

    function displayQueriesPanel(hcbo, eventStruct)
        qPanel = uipanel(f, 'Title', 'queries', 'position',[.01 .2 2 .8]);
        qscrollbar=uicontrol('parent', qPanel, 'style','slider','units','normalized','position',[0 0 1 .05],'callback',{@qpanel_scroll, qPanel});
        startPosition = .03; height = .9;
        for l = 1:length(queries)
            queryString(l) = uicontrol('parent', qPanel,'units', 'normalized', 'Style', 'Text', 'String', queries{l}, 'position', [startPosition height .5 .08], 'HorizontalAlignment', 'left');
            queryCheckBox(l) = uicontrol('parent', qPanel, 'units', 'normalized','Style', 'checkbox', 'String', l, 'callback', {@checkbox, l},'position', [(startPosition-.02),(height+.033),.02,.05]);
            height = height - .05;
            if height >= 300
                height = .9;
                startPosition = 1;
            end
        end
        function qpanel_scroll(src, evt, qPanel)
            set(qPanel,'position',[-get(src,'value') .2 2 .8])
        end
        function checkbox(hcbo, eventStruct, numQuery)
            checkedQueries{end+1} = queries{numQuery};
        end
    end
    function getMoreSpecs(hcbo, eventStruct)
        set(WHERE_button, 'visible', 'off');
        pos = get(WHERE_button, 'position');
        set(WHERE_button, 'string', 'selected');
        mod_pos = [0, 0, -30, 0];
        mod = uicontrol(f, 'style', 'text', 'String', 'WHERE', 'position', pos+mod_pos);
        opt_pos = [50, 0, 20, 20];
        set(secondary_action, 'position', pos+opt_pos);
        set(secondary_action, 'String', tot_columns);
        set(secondary_action, 'visible', 'on');
        set(another_equals, 'visible', 'on');
        eq_pos = [160, 0, -70, 0];
        set(another_equals, 'position', pos+eq_pos);
    end
    function AQ(hcbo, eventStruct)
        close(f);
        SqlGUI(queries);
    end
    function drop(hcbo, eventStruct)
        if ~isempty(checkedQueries)
            w = warndlg(checkedQueries{:}, 'Pressing OK will drop following queries:');
            waitfor(w);
            queryIDX = find(ismember(queries, checkedQueries));
            for i = 1:length(queryIDX)
                queries(queryIDX(i)) = [];
            end
            if isempty(queries)
                %close(conn);
                AQ();
            else
                set(findall(f, 'type', 'uicontrol'), 'visible', 'off');
                onArray = [Drop, Another, DB_send, Cancel, CloseAll];
                set(onArray, 'visible', 'on');
                height = 0;
                for l = 1:length(queries)
                    queryString(l) = uicontrol(f, 'Style', 'Text', 'String', queries{l}, 'position', [50, 336-height, 550, 40], 'HorizontalAlignment', 'left');
                    queryCheckBox(l) = uicontrol(f, 'Style', 'checkbox', 'String', l, 'callback', {@checkbox, l},'position', [20, 350-height, 30, 40]);
                    height = height + 40;
                end
            end
        else
            %close(conn);
            AQ();
        end
    end
    function done(hcbo, eventStruct)
        if ~isempty(checkedQueries)
            for i = 1:length(checkedQueries)
                if sscanf(checkedQueries{i}, 'DELETE')
                    done = questdlg(['Your ', i, 'th query will delete items from the table. Are you sure you''d like to delete?'], ...
                        'Please confirm', 'Yes, permanently delete', 'no, don''t include the query', 'no, don''t include the query');
                    waitfor(done);
                    if strcmp(done, 'no, don''t include the query')
                        checkedQueries(i) = [];
                    end
                end
                sent = exec(conn, checkedQueries(i));
                if sent.Message
                    w1 = warndlg(sent.Message, 'Error in query');
                    waitfor(w1);
                    keyboard;
                else
                    if sscanf(checkedQueries{i}, 'SELECT')
                        if ~isempty(secondary_actions)
                            export2wsdlg('Save value to workspace', secondary_actions, sent.Data, 'Save returned variables to workspace?');
                        end
                    end
                    queryIDX = ismember(checkedQueries(i), checkedQueries);
                    checkedQueries(queryIDX) = [];
                end
            end
        else
            warndlg('No query checkboxes selected');
        end
    end
    function final_values_string =  mapToString(cond_map, v_or_c)
        temp_values_string = 'WHERE';
        cond_keys = keys(cond_map);
%         date_input_string = '';
        %         if strcmp(values_keys, 'date') %FIXME BASED ON WHAT YOU STORE
        %             if values_map('propixx');
        %                 [values_map('start_date'), values_map('end_date')] = adjustDatesIfPropixxOnly(values_map('start_date'), values_map('end_date')];
        %             end
        %             date_input_string = selectDateInputQuerySubstring(values_map('start_date'), values_map('end_date'), values_map('sole_date'));
        %             rm_keys = {'end_date', 'start_date', 'sole_date', 'propixx'};
        %             remove(values_map, rm_keys);
        %         end
        for i = 1:size(cond_map, 1)
            v = values(cond_map, cond_keys(i));
            if v_or_c
                temp_values_string = sprintf('%s ''%s''=''%s'',', temp_values_string, cond_keys{i}, v{1});
            else
                temp_values_string = sprintf('%s AND ''%s''=''%s'',', temp_values_string, cond_keys{i}, v{1});
            end
        end
        if ~v_or_c
            temp_values_string = temp_values_string(4:end);
        end
%         if ~isempty(date_input_string)
%             final_values_string = strcat(date_input_string, ' ', temp_values_string(1:end-1));
%         else
            final_values_string = temp_values_string(1:end-1);
%         end
    end
    function cell_string = cellToString(cell)
        for i = 1:length(cell)
            cell{i} = ['''' cell{i} '''' ', '];
        end
        cell_string_temp = cell2mat(cell');
        cell_string = cell_string_temp(1:end-2);
    end
end