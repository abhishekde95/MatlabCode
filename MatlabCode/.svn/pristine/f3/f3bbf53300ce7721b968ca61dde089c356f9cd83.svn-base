%MAIN WORKHORSE FOR EXAMINING LMTF DATA. ANY OTHERS ARE DEPRECATED. 
%Fundamentally, this function takes a file list from
%the database (for a single subject, collects the raw data from the files, 
%organizes it, and plots the green funnel shapes and residuals, testing the
%residuals for any unusually high or low points. To do so, set model to -1. 
%Additionally, this function can be used to test/plot green funnel shapes
%and residuals for model fits to the data. To use this functionality, set
%"model" to the model number of your choice. setting the model to "0"
%prematurely exits the function. 

function [oddresids, data] = iterateAndPlotFiles_modularPlusDB(flist, model)
oddresids = {}; suspiciousfiles = {};
% put that data into the correct format, by iterating over the list of .nex files.
[data, sanityCheck] = getLMTFrawdata(flist);
RFXYData = unique(data(:,[5 6]),'rows');
single_filename = flist{end};
subjectID = single_filename(1); %Maybe change to take more than single letter subject IDs?
% At this point "data" has n rows and 7 columns. The n rows correspond to n
% different stimulus directions (L,M,TF combinations) that were tested. The
% four columns are (1) L-cone contrast, (2) M-cone contrast, (3) Temporal
% frequency, (4) out of gamut (1 or 0), (5) X-coordinate of
% eccentricity, (6) Y-coordinate of eccentricity, and (7) file identifier
% -----------------------------------------
% Now we're transitioning from collecting the data from the various files
% to fitting the data with various models
for i = 1:size(RFXYData,1)
    Lecc = data(:,5) == RFXYData(i,1) & data(:,6) == RFXYData(i,2);
    if model==-1
        [modelParams, curr_error] = LMTFfitting(data(Lecc,:),0);
    elseif model == 0
        disp('currently there is no functionality corresponding to model 0');
        return;
    else
        %change the first line to load the modelstruct from the local
        %machine.
        global_modelstruct = load('C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Zack\IsoSamp\private\data\LMTF.mat', subjectID);
        global_modelParams = getfield(global_modelstruct, subjectID);
        modelParams = LMTF_global_to_local_model(global_modelParams.model, RFXYData(i,1)/10, RFXYData(i,2)/10, model);
    end
    %extract isodetection data and populate cell array IGNORE FOR NOW
    %isoData(iterator,:) = horzcat(RFXYData, fpar); IGNORE FOR NOW
    splitTitle = {num2str(RFXYData(i,:))};
    
    % plot data in funnels and name that plot as the txt file's name
    %rawDataPlot(data(Lecc,:), modelParams, splitTitle, 0, 0); %just the data cloud, useful for figure making purposes only
    %rawDataPlot(data(Lecc,:), modelParams, splitTitle, 0, 1); %just the data cloud, zoomed to lowest 10tfs, for figures only
    rawDataPlot(data(Lecc,:), modelParams, splitTitle, 1, 0); %data cloud + fit THIS IS USEFUL FOR REGULAR DATA COLLECTION
    %rawDataPlot(data(Lecc,:), modelParams, splitTitle, 1, 1); %data cloud + fit, zoomed to lowest 10tfs, for figures only
    
    residsPlot(data(Lecc,:), modelParams, splitTitle); %plots the residuals to the data
    
    % Looking for individual data file with unusually large (or small)
    % residuals. GDLH 8/27/16
    predr = LMTF_thresh_from_model(data(Lecc,1),data(Lecc,2),data(Lecc,3), modelParams); % vector lengths in cc space
    r = sqrt(data(Lecc,1).^2+data(Lecc,2).^2);
    resids = r - predr;
    fractionalresid = resids./r;
    %plot resids along a line - IGNORE FOR NOW
    %residsLinePlot(splitTitle, flist(data(Lecc, 7)), fractionalresid);
    %determine if resids are funky
    try
        file_idx = data(Lecc,7);
    catch
        keyboard;
    end
    Loog = data(Lecc,4);
    [file_idx(~Loog) fractionalresid(~Loog)];
    mn_frac_idx = [];
    for j = unique(file_idx)'
        mn_frac_idx = [mn_frac_idx; mean(fractionalresid(file_idx==j & ~Loog))];
        if (mean(fractionalresid(file_idx==j & ~Loog)) < -.5)
            disp(['File with unusually low residuals found: ',char(flist(j))])
            suspiciousfiles(length(suspiciousfiles)+1) = flist(j);
            r = fractionalresid(file_idx==j & ~Loog);
            disp(r);
            oddresids{length(oddresids)+1} = {'low', flist(j), r};
        elseif (mean(fractionalresid(file_idx==j & ~Loog)) > .5)
            disp(['File with unusually high residuals found: ',char(flist(j))])
            suspiciousfiles(length(suspiciousfiles)+1) = flist(j);
            r = fractionalresid(file_idx==j & ~Loog);
            oddresids{length(oddresids)+1} = {'high', flist(j), r};
        end
    end
end
%plot elliptical slices of data - if several datasets are stitched
%together, this will plot all of them on one graph. PROBABLY NOT FUNCTIONAL
%binnedDataPlot(data, 1,1,1, splitTitle{1}); %plot 4, plot with dots 
end