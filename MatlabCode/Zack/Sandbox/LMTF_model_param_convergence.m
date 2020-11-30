textfiles_path = nexfilepath('nexfilelists', 'Greg', 'LMTF');
text_filenames = dir(fullfile(textfiles_path, '*LMTF.txt'));
subjectIDs = cellfun(@(x) {x(1)}, {text_filenames.name}); % the code below assumes all subject IDs are unique
text_files = strcat(textfiles_path, filesep, {text_filenames.name});
subject_counter = 0;
for text_file = text_files
    subject_counter = subject_counter + 1;
    flist = flatten(fnamesFromTxt2(text_file{1}));
    data = [];
    models = zeros(12+1,length(flist)); % tack theta on the end
    for k = 1:length(flist)
        stro = nex2stro(findfile(flist{k}));
        Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
        Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
        Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
        Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
        
        [~,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
        questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
        tfs = stro.trial(init_stim_trial_idxs,Ltf);
        
        % Out of gamut checking
        funds = reshape(stro.sum.exptParams.fundamentals, [], 3);
        spds = reshape(stro.sum.exptParams.mon_spd, [], 3);
        spds = SplineSpd(linspace(380, 780, size(spds, 1))', spds, linspace(380, 780, size(funds, 1))');
        M = funds'*spds;
        bkgndrgb = stro.sum.exptParams.bkgndrgb;
        [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
        questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
        data = [data; questmodes tfs ~in_gamut']; %#ok<AGROW> % Lcc Mcc TF OOG
        
        Loog = logical(data(:,end));
        
        x = data(:,1);
        y = data(:,2);
        tf = data(:,3);
        
        LB = [0 0 1 1 .001 .001];
        UB = [100 1 20 20 .03 .03];
        
        initparams = [40 .1 9 3 .005 .002]; % Good for Greg
        options = optimset('Algorithm', 'active-set', 'MaxFunEvals', 50000, 'MaxIter', 50000, 'TolFun', 10^-8);
        thetas = linspace(0, pi/2, 12);
        toterr = zeros(size(thetas));
        fpar = [initparams, initparams];
        fpars = zeros(length(fpar), length(thetas));
        for n = 1:length(thetas)
            rotmat = [cos(thetas(n)) sin(thetas(n)); -sin(thetas(n)) cos(thetas(n))];
            % Rotating data clockwise = rotating fit axis counterclockwise
            xytmp = [x y]*rotmat';
            [fpar,fv] = fmincon(@(params) tf_fiterr2(params, [xytmp(:,1) xytmp(:,2) tf], Loog), fpar, ...
                [], [], [], [], [LB LB], [UB UB], [], options);
            toterr(n) = fv;
            fpars(:,n) = fpar;
        end
        [~,bestrotidx] = min(toterr);
        models(:,k) = [fpars(:,bestrotidx); thetas(bestrotidx)];
    end
    s = struct(subjectIDs{subject_counter}, models);
    if subject_counter == 1
        save('LMTF_all_models.mat', '-struct', 's', subjectIDs{subject_counter});
    else
        save('LMTF_all_models.mat', '-struct', 's', subjectIDs{subject_counter}, '-append');
    end
end
