%Called after every lmtf file is collected to add each file to database
%online. SUPER important - do not delete this or related files!!
function DBPostPopGUI()
%check for file in local plexon D
origPath = 'D:\PlexonData';
%origPath = 'N:\NexFiles\Emily\Utu\LMTF'; %THIS PATH IS FOR TESTING PURPOSES ONLY
newPath = 'N:\NexFiles\Lisa';
fileGrabber(origPath, newPath, 1);
end