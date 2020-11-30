% GDLH modified 8/27/16 to pass back file identifier number in the last
% column of "data"
% EG modified 11/29/16 to just take flist, and to output a reasonable
% list of problem files plus information regarding the problem.
% 5/8/18 - fundamental part of iterateAndPlotFiles_modularPlusDB. 
function [data, problemFile, bkgndrgb, M] = getLMTFrawdata(flist)
problemFile = [];
for i = 1:length(flist)
    try
        stro = nex2stro(findfile(flist{i})); % <-- process the information in the nex file and put it into the "stro" structure
        % Below, just figuring out what information is in what column
        Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
        Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
        Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
        Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
        
        % Getting the threshold points
        [~,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last'); %check on this
        questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
        tfs = stro.trial(init_stim_trial_idxs,Ltf);
        
        %getting eccentricities
        stim_x = repmat(stro.sum.exptParams.stim_x, size(tfs,1), 1);
        stim_y = repmat(stro.sum.exptParams.stim_y, size(tfs,1), 1);
        
        % Out of gamut checking
        funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
        if (size(stro.sum.exptParams.mon_spd,1) == 303)
            spds = SplineSpd([380:4:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
        else
            spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
        end
        M = funds'*spds;
        bkgndrgb = stro.sum.exptParams.bkgndrgb;
        if sum(isnan(questmodes))>0
            problemFile = [problemFile; flist(i) questmodes bkgndrgb];
            disp('nan questmodes');
            continue;
        end
        [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
        questmodes(~in_gamut,:) = questmodes(~in_gamut,:).*repmat(scalar(~in_gamut)', 1, 2);
        %removing files where tfs under 30 are oog
        if tfs(find(~in_gamut)) <= 40
            problemFile = [problemFile; flist(i) tfs in_gamut];
            continue;
        else
            if exist('data', 'var')
                data = [data; questmodes tfs ~in_gamut' abs(stim_x) stim_y repmat(i,length(questmodes),1)];
            else
                data = [questmodes tfs ~in_gamut' abs(stim_x) stim_y repmat(i,length(questmodes),1)];  
            end
        end
    catch
        keyboard;
        problemFile = [problemFile; flist(i) questmodes repmat(i,length(questmodes),1)];
        continue;
    end
end
end