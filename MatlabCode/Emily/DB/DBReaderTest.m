%trying to read in data from different forms of DB via matlab. First few
%attempts are using the excel sheet method, the last attempt is using the
%SQL method. The fields are as follows: File ID, Paradigm name, RFX, RFY,
%Quality, Notes, Neuron #. Not really useful - can easily be deleted.
%% Section 1: Excel
filename = 'C:\Users\emily.gelfand\Desktop\ApolloDB_sample.csv';
%not successful way to read in csv data - only reads ints
%A = csvread(filename);
%disp(A);

%import from excel database - functional
[num,txt,raw] = xlsread(filename);
disp(raw); %this is the full table, including headers, that we're looking for. 
%to read by content without reading the whole sheet (i.e. "only where
%neuron = 1"), you have to create a function. 
[idk, irdk, raw, idx] = xlsread(filename, 'ApolloDB_sample', '','', @selectNeuron);


%trying to import another way, haven't played with enough
C = readtable(filename);
disp(C); % view whole table, including headers
%no way to use this to read by content type (i.e. "only where neuron = 1") without reading the whole sheet,
%but you can read in data ranges (i.e. "first two rows");
%% Section 2: MySQL

%import SQL database
conn = database('HorwitzTest','admin','adminVector123','Vendor','MySql','Server','localhost'); %connect to database
curs = exec(conn, 'SELECT * FROM Apollo WHERE paradigmID='); %view whole table
curs2 = fetch(curs);
fullDataTableApollo = curs2.Data;
disp(fullDataTableApollo);
close(curs); close(curs2); %close the cursor you're finished with - to avoid conflicts

nextcurs = exec(conn, 'SELECT * FROM Apollo WHERE neuron=1'); %view one neuron's data
nextcurs2 = fetch(nextcurs);
partDataTableApollo = nextcurs2.Data;
disp(partDataTableApollo);
close(nextcurs); close(nextcurs2); %same as closes from before

close(conn); %exit database connection