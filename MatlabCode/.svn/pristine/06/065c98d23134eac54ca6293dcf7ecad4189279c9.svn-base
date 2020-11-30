datapath = 'N:\NexFiles';
cd(datapath);
cdir = dir;
subdirs = {};
alldone = 0;
while ~alldone
	for c = 1:size(cdir,1)
		%filenames = ??
		%DBinsert(filenames);
		if (cdir(c).isdir && ~strncmp(cdir(c).name, '.',1)) %and we haven't seen it before
			%save the name of the directory
			%set cdir to the new directory name
			%cd into directory
			cdir = Outermost_dirstruct(o);
			subdirs{length(subdirs)+1} = (cdir.name);
			cd(cdir);
		else
			continue;
		end
	end
	if (strncmp(cdir.name,datapath)) %check this
		for b = 1:size(cdir,1)
			ispresent = find(cellfun(@(s) ~isempty(strfind(cdir(b).name, s)), subdirs));
			if ispresent
				continue;
			else
				cd(cdir(b));
				break; %this should break out of this for loop and go up to the top of the while loop, thus starting the first for loop over again.
			end
			alldone = 1;
		end
	else
		cd('../');
	end
end

Function DBinsert(filenames)
	%connect to db
	for f = 1:length(filenames)
		%check if it's nex
		%check if it's openable
		%check if it's useful
		DBtables(filename, conn)
	end
	%close db connection

Function DBtables(filename, conn)
		%if tablename not in DB
			%create table
			%send to DB
		%end
		filedumper(filename, conn)

Function filedumper(filename, conn)
	%if filename not in table
		%collect info
		%send to DB
	%else
		%store info
	%end