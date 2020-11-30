%%
%SMurray_p task



clear all



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose how many files to analyze
num_files = 3; %'1' to analyze one file
               %'2' to analyze multiple files for one day of data collection
               %'3' to analyze multiple days of data collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if num_files == 1
    stro = nex2stro
elseif num_files == 2
    [fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
    day = fname([1:7 11:14]);
    folder = pathname;
    allfiles = dir(folder);
    rowfiles = 1;
    for i = 1:size(allfiles,1)
        if size(allfiles(i).name,2) == 14
            if allfiles(i).name([1:7 11:14]) == day
                allblocks(rowfiles,:) = allfiles(i).name;
                rowfiles = rowfiles + 1;
            end
        end
    end
    row = 1;
    rowC = 1;
    for k = 1:size(allblocks,1)
        stro1 = nex2stro(strcat(folder, allblocks(k,:)));
        bkgnd = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'bkgnd'));
        if bkgnd == 2 %gray
            stro.sum = stro1.sum;
            stro.trial(row:row+length(stro1.trial)-1,:) = stro1.trial(1:length(stro1.trial),:);
            stro.ras(row:row+length(stro1.ras)-1,:) = stro1.ras(1:length(stro1.ras),[1:2 end]);
            row = row+length(stro1.trial);
        elseif bkgnd == 0 %corridor
            stroC.sum = stro1.sum;
            stroC.trial(rowC:rowC+length(stro1.trial)-1,:) = stro1.trial(1:length(stro1.trial),:);
            stroC.ras(rowC:rowC+length(stro1.ras)-1,:) = stro1.ras(1:length(stro1.ras),[1:2 end]);
            rowC = rowC+length(stro1.trial);
        end
    end
elseif num_files == 3
    [fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
    folder = pathname;
    allfiles = dir(folder);
    rowfiles = 1;
    for i = 1:size(allfiles,1)
        if size(allfiles(i).name,2) == 14
            if allfiles(i).name(11:14) == fname(11:14)
                allblocks(rowfiles,:) = allfiles(i).name;
                rowfiles = rowfiles + 1;
            end
        end
    end
    daynames = unique(allblocks(:, 1:7), 'rows');
    row = 1;
    rowC = 1;
    rowday = ones(size(daynames,1), 1);
    rowdayC = ones(size(daynames,1), 1);
    for k = 1:size(allblocks,1)
        stro1 = nex2stro(strcat(folder, allblocks(k,:)));
        bkgnd = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'bkgnd'));
        if bkgnd == 2 %gray
            %combining all trials with this background
            stro.sum = stro1.sum;
            stro.trial(row:row+length(stro1.trial)-1,:) = stro1.trial(1:length(stro1.trial),:);
            stro.ras(row:row+length(stro1.ras)-1,:) = stro1.ras(1:length(stro1.ras),[1:2 end]);
            row = row+length(stro1.trial);
            %separating data into separate days of data collection
            clear sacc_time_nearTEMP stim_or_nearTEMP stim_or_farTEMP
            sacc_time_nearTEMP = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'sacc_time_near'));
            stim_or_nearTEMP = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'stim_or_near'));
            stim_or_farTEMP = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'stim_or_far'));
            whichday = strcmp(allblocks(k,1:7), cellstr(daynames));
            dn = find(whichday);
            stroDay(dn).sacc_time_nearDAY(rowday(dn):rowday(dn)+length(sacc_time_nearTEMP)-1,1) = sacc_time_nearTEMP;
            stroDay(dn).stim_or_nearDAY(rowday(dn):rowday(dn)+length(stim_or_nearTEMP)-1,1) = stim_or_nearTEMP;
            stroDay(dn).stim_or_farDAY(rowday(dn):rowday(dn)+length(stim_or_farTEMP)-1,1) = stim_or_farTEMP;
            rowday(dn) = rowday(dn)+length(sacc_time_nearTEMP);
        elseif bkgnd == 0 %corridor
            %combining all trials with this background
            stroC.sum = stro1.sum;
            stroC.trial(rowC:rowC+length(stro1.trial)-1,:) = stro1.trial(1:length(stro1.trial),:);
            stroC.ras(rowC:rowC+length(stro1.ras)-1,:) = stro1.ras(1:length(stro1.ras),[1:2 end]); 
            rowC = rowC+length(stro1.trial);
            %separating data into separate days of data collection
            clear sacc_time_nearTEMP stim_or_nearTEMP stim_or_farTEMP
            sacc_time_nearTEMP = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'sacc_time_near'));
            stim_or_nearTEMP = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'stim_or_near'));
            stim_or_farTEMP = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'stim_or_far'));
            whichday = strcmp(allblocks(k,1:7), cellstr(daynames));
            dn = find(whichday);
            stroDayC(dn).sacc_time_nearDAY(rowdayC(dn):rowdayC(dn)+length(sacc_time_nearTEMP)-1,1) = sacc_time_nearTEMP;
            stroDayC(dn).stim_or_nearDAY(rowdayC(dn):rowdayC(dn)+length(stim_or_nearTEMP)-1,1) = stim_or_nearTEMP;
            stroDayC(dn).stim_or_farDAY(rowdayC(dn):rowdayC(dn)+length(stim_or_farTEMP)-1,1) = stim_or_farTEMP;
            rowdayC(dn) = rowdayC(dn)+length(sacc_time_nearTEMP); 
        end
    end
end


%%
%task and stimulus parameters

bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
ntrials = size(stro.trial,1);
starttimes = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
sacc_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sacc_time'));
sacc_time_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sacc_time_near'));
sacc_time_far = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sacc_time_far'));
stim_or_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or_near'));
stim_or_far = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or_far'));
correct_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct_resp'));

uniquenf = unique([stim_or_near stim_or_far],'rows');
uniquen = unique(stim_or_near,'rows');
uniquef = unique(stim_or_far,'rows');

if exist('stroC','var')
    ntrialsC = size(stroC.trial,1);
    starttimesC = stroC.ras(:,strcmp(stroC.sum.rasterCells,'anlgStartTime'));
    stimon_tC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'stim_on'));
    stimoff_tC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'stim_off'));
    targon_tC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'targ_on'));
    sacc_tC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'sacc_time'));
    sacc_time_nearC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'sacc_time_near'));
    sacc_time_farC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'sacc_time_far'));
    stim_or_nearC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'stim_or_near'));
    stim_or_farC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'stim_or_far'));
    correct_tC = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'correct_resp'));
    
    uniquenfC = unique([stim_or_nearC stim_or_farC],'rows');
    uniquenC = unique(stim_or_nearC,'rows');
    uniquefC = unique(stim_or_farC,'rows');
end
    

%%
%reaction time histogram
if all(isnan(sacc_t))
    saccnear = find(~isnan(sacc_time_near));
    sacc_t(saccnear) = sacc_time_near(saccnear);
    saccfar = find(~isnan(sacc_time_far));
    sacc_t(saccfar) = sacc_time_far(saccfar);
end
rt = sacc_t - targon_t;
rt = rt*1000; %in milliseconds

if exist('stroC','var')
    if all(isnan(sacc_tC))
        saccnearC = find(~isnan(sacc_time_nearC));
        sacc_tC(saccnearC) = sacc_time_nearC(saccnearC);
        saccfarC = find(~isnan(sacc_time_farC));
        sacc_tC(saccfarC) = sacc_time_farC(saccfarC);
    end
    rtC = sacc_tC - targon_tC;
    rtC = rtC*1000; %in milliseconds
end

figure
hold on
[freqrt, binrt] = hist(rt);
if bkgnd == 2
    bar(binrt, freqrt, 1, 'k')
elseif bkgnd == 0
    bar(binrt, freqrt, 1, 'r')
end
if exist('stroC','var')
    [freqrtC, binrtC] = hist(rtC);
    bar(binrtC, freqrtC, 1, 'r')
    legend('gray', 'corridor')
end
xlabel('reaction time (ms)')
ylabel('trials')
hold off


%%
%average reaction time per inter-stimulus size difference

stim_diff = abs(stim_or_near - stim_or_far); 
stim_diff = round(stim_diff*10) / 10; %rounded to hundredth of a degree
uniquediff = unique(stim_diff); %unique inter-stimulus size differences
stim_diff(:,2) = rt;
rt_diff = NaN(size(uniquediff,1), 3);
for i = 1:size(uniquediff, 1)
    each = uniquediff(i);
    rt_diff(i, 1) = each;
    eachstim = stim_diff(:,1) == each;
    allrt = stim_diff(eachstim,2);
    rt_diff(i, 2) = mean(allrt);
    rt_diff(i, 3) = std(allrt) / sqrt(length(allrt)); %SEM
end
rt_diff(:,1) = rt_diff(:,1) ./ 10; %degrees

if exist('stroC','var')
    stim_diffC = abs(stim_or_nearC - stim_or_farC);
    stim_diffC = round(stim_diffC*10) / 10; %rounded to hundredth of a degree
    uniquediffC = unique(stim_diffC); %unique inter-stimulus size differences
    stim_diffC(:,2) = rtC;
    rt_diffC = NaN(size(uniquediffC,1), 3);
    for i = 1:size(uniquediffC, 1)
        eachC = uniquediffC(i);
        rt_diffC(i, 1) = eachC;
        eachstimC = stim_diffC(:,1) == eachC;
        allrtC = stim_diffC(eachstimC,2);
        rt_diffC(i, 2) = mean(allrtC);
        rt_diffC(i, 3) = std(allrtC) / sqrt(length(allrtC)); %SEM
    end
    rt_diffC(:,1) = rt_diffC(:,1) ./ 10; %degrees  
end

figure
hold on
if bkgnd == 2
    errorbar(rt_diff(:,1), rt_diff(:,2), rt_diff(:,3), 'k')
elseif bkgnd == 0
    errorbar(rt_diff(:,1), rt_diff(:,2), rt_diff(:,3), 'r')
end
if exist('stroC','var')
    errorbar(rt_diffC(:,1), rt_diffC(:,2), rt_diffC(:,3), 'r')
    legend('gray', 'corridor')
end
xlabel('inter-stim size difference (degrees)')
ylabel('average reaction time (ms)')
hold off


%%
%GRAY BACKGROUND ONLY
%overall performance based on inter-stimulus size difference

bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
if bkgnd == 2
    stim_diff(:,3) = correct_t;
    perf_diff = NaN(size(uniquediff,1), 2);
    for i = 1:size(uniquediff, 1)
        each = uniquediff(i);
        perf_diff(i, 1) = each;
        eachstim = stim_diff(:,1) == each;
        alldiff = stim_diff(eachstim,3);
        allwrong = find(isnan(alldiff));
        perf_diff(i, 2) = ((length(alldiff) - length(allwrong)) / length(alldiff)) * 100; %percent correct
    end
    perf_diff(:,1) = perf_diff(:,1) ./ 10; %degrees
    
    figure
    hold on
    plot(perf_diff(:,1), perf_diff(:,2), '.-k')
    legend('gray', 'Location', 'East')
    xlabel('inter-stim size difference (degrees)')
    ylabel('% correct performance')
    hold off
end

%%
%performance based on stimulus location (near vs. far)

nfdiff = stim_or_near - stim_or_far; %how much bigger the near stim is than the far stim
nfdiff = round(nfdiff*10) / 10; %rounded to hundredth of a degree
diffunique = unique(nfdiff);
%determines if the saccade was made to "near"
if all(isnan(sacc_time_near))
    for i = 1:size(nfdiff,1)
        if nfdiff(i, 1) > 0 && correct_t(i,1) > 0
            nfdiff(i, 2) = 1; %near choice
        elseif nfdiff(i, 1) > 0 && isnan(correct_t(i,1))
            nfdiff(i, 2) = 0; %far choice
        elseif nfdiff(i, 1) < 0 && isnan(correct_t(i,1))
            nfdiff(i, 2) = 1; %near choice
        elseif nfdiff(i, 1) < 0 && correct_t(i,1) > 0
            nfdiff(i, 2) = 0; %far choice
        end
    end
else
    nfdiff(:,2) = 0;
    saccnear = find(~isnan(sacc_time_near));
    nfdiff(saccnear,2) = 1;
end
nearchoice = nfdiff(:,2) == 1;
for i = 1:size(diffunique,1)
    d = diffunique(i,1) == nfdiff(:,1);
    total = length(find(d));
    dn = d & nearchoice;
    per = length(find(dn));
    diffunique(i,2) = per/total * 100;
end
diffunique(:,1) = diffunique(:,1) ./ 10; %degrees




figure
hold on
if bkgnd == 2
    plot(diffunique(:,1), diffunique(:,2), '.-k')
elseif bkgnd == 0
    plot(diffunique(:,1), diffunique(:,2), '.-r')
end
xlabel('near radius > far radius (degrees)')
ylabel('% near chosen')
%percentage 'near' stimulus was chosen across all trials
nearchosen = 0;
for i = 1:size(sacc_time_near,1)
    if ~isnan(sacc_time_near(i,1))  %saccade to near
        nearchosen = nearchosen + 1;
    end
end
percentnearchosen = nearchosen / size(sacc_time_near,1) * 100;

if exist('stroC','var')
    nfdiffC = stim_or_nearC - stim_or_farC; %how much bigger the near stim is than the far stim
    nfdiffC = round(nfdiffC*10) / 10; %rounded to hundredth of a degree
    diffuniqueC = unique(nfdiffC);
    %determines if the saccade was made to "near"
    nfdiffC(:,2) = 0;
    saccnearC = find(~isnan(sacc_time_nearC));
    nfdiffC(saccnearC,2) = 1;
    nearchoiceC = nfdiffC(:,2) == 1;
    for i = 1:size(diffuniqueC,1)
        dC = diffuniqueC(i,1) == nfdiffC(:,1);
        totalC = length(find(dC));
        dnC = dC & nearchoiceC;
        perC = length(find(dnC));
        diffuniqueC(i,2) = perC/totalC * 100;
    end
    diffuniqueC(:,1) = diffuniqueC(:,1) ./ 10; %degrees
    plot(diffuniqueC(:,1), diffuniqueC(:,2), '.-r')
    %CORRIDOR BACKGROUND: percentage 'near' stimulus was chosen across all trials
    nearchosenC = 0;
    for i = 1:size(sacc_time_nearC,1)
        if ~isnan(sacc_time_nearC(i,1))  %saccade to near
            nearchosenC = nearchosenC + 1;
        end
    end
    percentnearchosenC = nearchosenC / size(sacc_time_nearC,1) * 100;
end

if exist('stroC','var')
    legend(['gray: near chosen ' num2str(round(percentnearchosen)) '%'], ... 
           ['corridor: near chosen ' num2str(round(percentnearchosenC)) '%'], 'Location', 'NorthWest')
else
    legend(['near chosen ' num2str(round(percentnearchosen)) '%'], 'Location', 'NorthWest')
end
plot([0 0], [0 100], '--k')
plot([min(diffunique(:,1)) max(diffunique(:,1))], [50 50], '--k')
axis([min(diffunique(:,1)) max(diffunique(:,1)) 0 100]);
hold off


%%
%if num_files == 3
%plot of performance based on stimulus location (near vs. far)
%but taking the average across days (n = # days) instead of trials (n = # trials)


if num_files == 3
    row = 1; %assuming size of 'diffuniqueTmp' and 'diffuniqueTmpC' equal per day of data
    percentnearchosenDays = NaN(size(stroDay,2), 1);
    percentnearchosenDaysC = NaN(size(stroDayC,2), 1);
    for d = 1:size(stroDay,2)
        clear nfdiffTmp nfdiffTmpC diffuniqueTmp diffuniqueTmpC saccnearTmp saccnearTmpC 
        clear nearchoiceTmp nearchoiceTmpC
        nfdiffTmp = stroDay(d).stim_or_nearDAY - stroDay(d).stim_or_farDAY; %how much bigger the near stim is than the far stim
        nfdiffTmpC = stroDayC(d).stim_or_nearDAY - stroDayC(d).stim_or_farDAY;
        nfdiffTmp = round(nfdiffTmp*10) / 10; %rounded to hundredth of a degree
        nfdiffTmpC = round(nfdiffTmpC*10) / 10; 
        diffuniqueTmp = unique(nfdiffTmp); %unique size differences between the two stimuli
        diffuniqueTmpC = unique(nfdiffTmpC);
        %determines if the saccade was made to "near"
        nfdiffTmp(:,2) = 0;
        nfdiffTmpC(:,2) = 0;
        saccnearTmp = find(~isnan(stroDay(d).sacc_time_nearDAY));
        saccnearTmpC = find(~isnan(stroDayC(d).sacc_time_nearDAY));
        nfdiffTmp(saccnearTmp,2) = 1;
        nfdiffTmpC(saccnearTmpC,2) = 1;
        nearchoiceTmp = nfdiffTmp(:,2) == 1;
        nearchoiceTmpC = nfdiffTmpC(:,2) == 1;
        for i = 1:size(diffuniqueTmp,1) %assuming size of 'diffuniqueTmpC' is equal per day of data
            clear dTmp dTmpC dnTmp dnTmpC
            dTmp = diffuniqueTmp(i,1) == nfdiffTmp(:,1);
            dTmpC = diffuniqueTmpC(i,1) == nfdiffTmpC(:,1);
            totalTmp = length(find(dTmp));
            totalTmpC = length(find(dTmpC));
            dnTmp = dTmp & nearchoiceTmp;
            dnTmpC = dTmpC & nearchoiceTmpC;
            perTmp = length(find(dnTmp));
            perTmpC = length(find(dnTmpC));
            diffuniqueTmp(i,2) = perTmp/totalTmp * 100;
            diffuniqueTmpC(i,2) = perTmpC/totalTmpC * 100;
        end
        diffuniqueTmp(:,1) = diffuniqueTmp(:,1) ./ 10; %degrees
        diffuniqueTmpC(:,1) = diffuniqueTmpC(:,1) ./ 10; %degrees
        %percentage 'near' stimulus was chosen across all trials
        %for gray background first
        nearchosenTmp = 0;
        for i = 1:size(stroDay(d).sacc_time_nearDAY,1)
            if ~isnan(stroDay(d).sacc_time_nearDAY(i,1))  %saccade to near
                nearchosenTmp = nearchosenTmp + 1;
            end
        end
        percentnearchosenTmp = nearchosenTmp / size(stroDay(d).sacc_time_nearDAY,1) * 100;
        %for corridor background next
        nearchosenTmpC = 0;
        for i = 1:size(stroDayC(d).sacc_time_nearDAY,1)
            if ~isnan(stroDayC(d).sacc_time_nearDAY(i,1))  %saccade to near
                nearchosenTmpC = nearchosenTmpC + 1;
            end
        end
        percentnearchosenTmpC = nearchosenTmpC / size(stroDayC(d).sacc_time_nearDAY,1) * 100;
        %store each day's data
        diffuniqueDays(row:row+size(diffuniqueTmp,1)-1,:) = diffuniqueTmp;
        diffuniqueDaysC(row:row+size(diffuniqueTmpC,1)-1,:) = diffuniqueTmpC;
        row = row + size(diffuniqueTmp,1);
        %store the percentage 'near' stimulus was chosen across all trials
        percentnearchosenDays(d,1) = percentnearchosenTmp;
        percentnearchosenDaysC(d,1) = percentnearchosenTmpC;
    end
    
    %find average of data and SEM across days of data (n = # days)
    x = unique(diffuniqueDays(:,1)); %assuming the same stimulus sizes used for both bkgnds
    avgDays = NaN(length(x), 3);
    avgDaysC = NaN(length(x), 3);
    for i = 1:length(x)
        dx = x(i,1) == diffuniqueDays(:,1);
        dxC = x(i,1) == diffuniqueDaysC(:,1);
        totalx = length(find(dx));
        totalxC = length(find(dxC));
        %gray background first
        avgDays(i,1) = x(i);
        avgDays(i,2) = mean(diffuniqueDays(dx,2));
        avgDays(i,3) = std(diffuniqueDays(dx,2)) / sqrt(totalx); 
        %corridor background next
        avgDaysC(i,1) = x(i);
        avgDaysC(i,2) = mean(diffuniqueDaysC(dxC,2));
        avgDaysC(i,3) = std(diffuniqueDaysC(dxC,2)) / sqrt(totalxC); 
    end
    figure;
    hold on;
    errorbar(avgDays(:,1), avgDays(:,2), avgDays(:,3), '.-k')
    errorbar(avgDaysC(:,1), avgDaysC(:,2), avgDaysC(:,3), '.-r')
    legend(['gray: near chosen ' num2str(round(mean(percentnearchosenDays))) '%'], ...
        ['corridor: near chosen ' num2str(round(mean(percentnearchosenDaysC))) '%'], 'Location', 'NorthWest')
    plot([0 0], [0 100], '--k')
    plot([min(x) max(x)], [50 50], '--k')
    axis([min(x) max(x) 0 100]);
    xlabel('near radius > far radius (degrees)')
    ylabel('% near chosen')
    hold off
end



%%
%display percentage 'near' stimulus was randomly assigned to be the target when stimuli were equal size
%displays for GRAY BACKGROUND ONLY

nearequal = 0;
for i = 1:size(nfdiff,1)
    if     nfdiff(i,1) == 0 && ~isnan(sacc_time_near(i,1)) && ~isnan(correct_t(i,1)) %saccade to near, correct
        nearequal = nearequal + 1;
    elseif nfdiff(i,1) == 0 &&  isnan(sacc_time_near(i,1)) &&  isnan(correct_t(i,1)) %saccade to far, incorrect
        nearequal = nearequal + 1;
    end
end
equal = find(nfdiff(:,1) == 0);
equal = length(equal);
percentnear = nearequal / equal * 100;
if bkgnd == 2
    disp(['GRAY: percent near stim assigned as target when stim equal size: ' num2str(percentnear)]);
end


%%
%display percentage 'near' stimulus was chosen when stimuli were equal size

nearchosenequal = 0;
for i = 1:size(nfdiff,1)
    if nfdiff(i,1) == 0 && ~isnan(sacc_time_near(i,1))  %saccade to near
        nearchosenequal = nearchosenequal + 1;
    end
end
percentnearchosenequal = nearchosenequal / equal * 100;
if bkgnd == 2
    disp(['GRAY: percent near stim chosen when stim equal size: ' num2str(percentnearchosenequal)]);
elseif bkgnd == 0
    disp(['CORRIDOR: percent near stim chosen when stim equal size: ' num2str(percentnearchosenequal)]);
end

if exist('stroC','var')
    nearchosenequalC = 0;
    for i = 1:size(nfdiffC,1)
        if nfdiffC(i,1) == 0 && ~isnan(sacc_time_nearC(i,1))  %saccade to near
            nearchosenequalC = nearchosenequalC + 1;
        end
    end
    equalC = find(nfdiffC(:,1) == 0);
    equalC = length(equalC);
    percentnearchosenequalC = nearchosenequalC / equalC * 100;
    disp(['CORRIDOR: percent near stim chosen when stim equal size: ' num2str(percentnearchosenequalC)]);
    
end


%%
%performance in each combination of stimulus sizes

eachcombo(:,1) = round(stim_or_near*10) ./ 10; %rounded to hundredth of a degree
eachcombo(:,2) = round(stim_or_far*10) ./ 10;
eachcombo = eachcombo ./ 10; %degrees
uniquenear = unique(eachcombo(:,1));
uniquefar = unique(eachcombo(:,2));
z = NaN(length(uniquenear), length(uniquefar));
for i = 1:length(uniquenear)
    for j = 1:length(uniquefar)
        nearmatch = eachcombo(:,1) == uniquenear(i);
        farmatch = eachcombo(:,2) == uniquefar(j);
        bothmatch = nearmatch & farmatch;
        totalall = nfdiff(bothmatch, 2);
        totalnear = find(totalall == 1);
        z(j, i) = length(totalnear) / length(totalall) *100;
    end
end

figure
if bkgnd == 2
    title('GRAY')
elseif bkgnd == 0
    title('CORRIDOR')
end
hold on
grid on
surf(uniquenear, uniquefar, z)
axis([min(uniquenear) max(uniquenear) min(uniquefar) max(uniquefar) 0 100])
xlabel('near radius')
ylabel('far radius')
colorbar
hold off


if exist('stroC','var')
    eachcomboC(:,1) = round(stim_or_nearC*10) ./ 10; %rounded to hundredth of a degree
    eachcomboC(:,2) = round(stim_or_farC*10) ./ 10;
    eachcomboC = eachcomboC ./ 10; %degrees
    uniquenearC = unique(eachcomboC(:,1));
    uniquefarC = unique(eachcomboC(:,2));
    zC = NaN(length(uniquenearC), length(uniquefarC));
    for i = 1:length(uniquenearC)
        for j = 1:length(uniquefarC)
            nearmatchC = eachcomboC(:,1) == uniquenearC(i);
            farmatchC = eachcomboC(:,2) == uniquefarC(j);
            bothmatchC = nearmatchC & farmatchC;
            totalallC = nfdiffC(bothmatchC, 2);
            totalnearC = find(totalallC == 1);
            zC(j, i) = length(totalnearC) / length(totalallC) *100;
        end
    end
    
    figure
    title('CORRIDOR')
    hold on
    grid on
    surf(uniquenearC, uniquefarC, zC)
    axis([min(uniquenearC) max(uniquenearC) min(uniquefarC) max(uniquefarC) 0 100])
    xlabel('near radius')
    ylabel('far radius')
    colorbar
    hold off
end