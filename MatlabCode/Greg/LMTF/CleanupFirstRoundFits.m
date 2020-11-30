% Fitting each dataset (data at each retinal location) using as an initial
% guess the best fit from each other data sets. Continue until none of the
% fits improve. From GrantBrainStorming Section 4.16
%
% called by both LMTF_generate_module_data.m and LMTFpop.m (section 5.1)
    function [outmodels, fvs] = CleanupFirstRoundFits(uniqueXYs, data, inmodels, fvs)
        options = optimset('Algorithm', 'sqp', 'MaxFunEvals', 5e4,'MaxIter', 5e4, 'TolFun', 1e-8, 'Display', 'off');
        paramList = MakeParamList(boundsSorter('frf'));
        LB = paramList(1:2:length(paramList));
        UB = paramList(2:2:length(paramList));
        initialguessidxs = 1:size(uniqueXYs,1);
        keepingtrack = 1; % Just need to start with some value in keeping track, it's cleared two lines below anyway
        % keepingtrack keeps track of which set of model parameters, used as an
        % initial guess, improves which model fit.
        waitbar_h = waitbar(0,'Please wait...');
        while ~isempty(keepingtrack)
            keepingtrack = [];
            for ii = 1:size(uniqueXYs,1) % data comes from retinal location i
                for jj = initialguessidxs % initial guess comes from model fit at location j
                    fractionalwaythrough = ((ii-1)*length(initialguessidxs)+find(jj==initialguessidxs))/(size(uniqueXYs,1)*length(initialguessidxs));
                    waitbar(fractionalwaythrough,waitbar_h);
                    Lecc = all(data(:,[5 6]) == repmat(uniqueXYs(ii,:),size(data,1),1),2);
                    % Just pulling out the data that's at the right retinal position (ecc).
                    try
                    [model,fv] = fmincon(@(params) tf_fiterr2(params,  data(Lecc,1:3),data(Lecc,4)), inmodels(:,jj),...
                        [],[],[],[],LB,UB,[],options);
                    catch
                        keyboard
                    end
                    if (fv < fvs(ii))
                        disp('--------------------------------');
                        fprintf('Fit to data at (%d, %d) is improved by guessing model params from (%d, %d)\n',uniqueXYs(ii,:),uniqueXYs(jj,:));
                        fprintf('Fitting error decreased from %d to %d\n',fvs(ii),fv)
                        disp('--------------------------------');
                        inmodels(:,ii) = model;
                        fvs(ii) = fv;
                        keepingtrack = [keepingtrack; ii jj];
                    end
                end
            end
            % The first column of keepingtrack contains models that were updated in
            % the last round. Use these as initial guesses in nest round.
            if (~isempty(keepingtrack))
                initialguessidxs = unique(keepingtrack(:,1))';
                disp('Improved model fits at retinal locations');
                disp(uniqueXYs(initialguessidxs,:));
            end % otherwise, we're done
        end
        outmodels = inmodels;
        close(waitbar_h);
    end
    
    function paramList = boundsSorter(model_mode)
        bounds.xi = [0 500];
        bounds.zeta = [0 1];
        bounds.n = [1 40];
        bounds.delta_n = [0 10];
        bounds.logtau = [-3 -1];
        bounds.logkappa = [log10(1.00001) log10(2)];
        bounds.theta = [0 pi/2];
        bounds.a = [-10 10];
        base_bounds_early = [bounds.xi bounds.zeta bounds.n bounds.delta_n bounds.logtau bounds.logkappa bounds.theta];
        base_bounds_later = [base_bounds_early(2:end-1) base_bounds_early(2:end)];
        switch model_mode
            case 'frf'
                paramList = {[base_bounds_early(1:end-2) base_bounds_early]};
            case 'mode1d'
                paramList = {base_bounds_early(1:end-1)};
            case 'mode1'
                paramList{1} = base_bounds_later(1:end-1);
                paramList{2} = repmat([bounds.xi; bounds.xi; bounds.theta],1,nuniqueXYs);
            case 'mode1p1'
                paramList{1} = [base_bounds_early(2,3:end-1), base_bounds_early(2:end)];
                paramList{2} = repmat([bounds.xi; bounds.xi; bounds.n],1,nuniqueXYs);
            case 'mode0'
                paramList{1} = base_bounds_later;
                paramList{2} = repmat([bounds.xi; bounds.xi],1,nuniqueXYs);
            case 'mode2'
                paramList{1} = base_bounds_later;
                paramList{2} = repmat([bounds.a],6,nuniqueXYs);
            case 'mode4'
                paramList{1} = base_bounds_later;
                paramList{2} = repmat([bounds.a],8,nuniqueXYs);
            otherwise %mode3, mode5
                paramList{1} = base_bounds_later;
                paramList{2} = repmat([bounds.a],7,nuniqueXYs);
        end
    end
    
        function list = MakeParamList(params)
        if size(params,2) == 2
            fixedparameters = params{1};
            variableparameters = params{2};
        else
            list = params{:};
            return;
        end
        if (size(fixedparameters,2) > 1) % forcing a column
            fixedparameters = fixedparameters';
        end
        suffix = [];
        for ii = 1:size(variableparameters,2)
            for jj = 1:size(variableparameters,1)
                suffix = [suffix; variableparameters(jj,ii)];
            end
        end
        list = [fixedparameters; suffix];
    end
