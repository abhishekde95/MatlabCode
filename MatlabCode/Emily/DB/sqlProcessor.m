%not useful - was originally planned as part of SqlFunctionGui, which was
%never completed as no longer needed. 
function sqlProcessor(DB_choice, table_name, values_map, conditions_map, type)
%TODO: Adjust to work with GUI input, confirm that nothing that gets passed
%to the server isn't in string form
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
tables = fetch(conn, 'SHOW FULL TABLES');
isTable = strcmp(table_name, tables{1,:});
if ~isTable
	disp('error, table_name invalid');
    return;
end
switch type
	case 1 
		sql_query = sprintf('UPDATE %s SET %s WHERE %s', table_name, values_string, conditions_string);
	case 2
		sql_query = sprintf('INSERT INTO %s IF NOT EXISTS %s', table_name, conditions_string);
	case 3
		sql_query = sprintf('DELETE FROM %s WHERE %s', table_name, conditions_string);
	case 4
		sql_query = sprintf('SELECT %s FROM %s WHERE %s', values_string, table_name, conditions_string);
	otherwise
		disp('error, not a known query type');
		return;
end
sendToDB = exec(conn, sql_query);
if sendToDB.Message
	keyboard;
end
close(conn);
end
%% DATES:
%%PROPIXX ONLY
function [start_date_adjusted, end_date_adjusted] =  adjustDatesIfPropixxOnly(start_date, end_date)
m_propixx_date = datenum(2016, 05, 27);
sql_propixx_date = '2016-05-27';
if strcmp(start_date, 'NULL') || convertDateToMatlab(start_date) < m_propixx_date
	newsd = sql_propixx_date;
else
	newsd = start_date;
end
if strcmp(end_date, 'NULL') || convertDateToMatlab(end_date) >= m_propixx_date
	newed = end_date;
else
	disp('error, end date is earlier than propixx start date');
	return;
end
[start_date_adjusted, end_date_adjusted] = [newsd, newed];
end
%%DEAL WITH DATES
function date_input_string =  selectDateInputQuerySubstring(start_date, end_date, sole_date)
no_start = strcmp(start_date, 'NULL');
no_end = strcmp(end_date, 'NULL');
if no_start
	if no_end
		date_input_string = 'NULL';
	elseif ~sole_date
		date_input_string = sprintf('AND recDate < %s', end_date);
	else
		date_input_string = sprintf('AND recDate = %s', end_date);
	end
elseif ~no_start && ~sole_date
	if no_end
		date_input_string = sprintf('AND recDate > %s', start_date);
	elseif convertDateToMatlab(start_date) < convertDateToMatlab(end_date)
		disp('error, end date is earlier than start date');
		return;
	else
		date_input_string = sprintf('AND recDate BETWEEN %s AND %s', start_date, end_date);
	end
else
	if no_end
		date_input_string = sprintf('AND recDate = %s', start_date);
	else
		date_input_string = sprintf('AND recDate IN(%s, %s)', start_date, end_date);
	end
end
end
function date =  convertDateToMatlab(init_date)
	date = datenum(init_date(1-4), init_date(6-7), init_date(9-10));
end