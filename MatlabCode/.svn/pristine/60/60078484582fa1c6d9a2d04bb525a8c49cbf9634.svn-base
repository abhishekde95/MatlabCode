%use iterateandplotfiles_modularplusdb instead. This is ancient.
subj = 'A'; x = 50; y = 0; startpath = 'C:/NO BACKUP/NexFiles';
filenames = pullFilenamesFromDB('lmtf', subj, x, y);
%rfxyDataTemp = []; 
data = [];
data = getLMTFrawdata(data, filenames, startpath);
%RFXYData = transpose(unique(rfxyDataTemp));
% At this point "data" has n rows and 6 columns. The n rows correspond to n
% different stimulus directions (L,M,TF combinations) that were tested. The
% four columns are (1) L-cone contrast, (2) M-cone contrast, (3) Temporal
% frequency, (4) out of gamut (1 or 0), (5) X-coordinate of
% eccentricity, and (6) Y-coordinate of eccentricity
% -----------------------------------------
% Now we're transitioning from collecting the data from the various files
% to fitting the data with various models
%modelParams = LMTFfitting(data);
modelParams = quickLMTFmodelfit(data(:,1), data(:,2), data(:,3), data(:,4));
splitTitle = 'fsdfsd'; %change later
%extract isodetection data and populate cell array
%isoData(iterator,:) = horzcat(RFXYData, fpar); IGNORE FOR NOW
rawDataPlot(data, modelParams, splitTitle);
%binnedDataPlot(data, splitTitle);
%plot elliptical slices of data - if several datasets are stitched
%together, this will plot all of them on one graph.
%CrossSections(data, besttheta);
% REGRESSION DIAGNOSTICS
%residsPlot(data, modelParams, splitTitle);
