function Analyze_data_KM(gl)
 clearvars;
% This is a new analysis file which is meant for analysing data for low sf (0 cycles/deg) orange, cyan, lime and magenta stimuli. 
% We want to see if we can replicate the results of Sakurai and Mullen, 2006 
NCONTRASTS = 50; % 50 for Rig 1 Data and Propixx
tot_corr_trials = zeros(NCONTRASTS,4);
tot_incorr_trials = zeros(NCONTRASTS,4);
if nargin == 0
    [fname, pathname] = uigetfile('*.mat', 'Select a data file','MultiSelect','on');
    if isequal(fname,0) || isequal(pathname,0)
        return
    end
    char_fname = char(fname);
    N = size(char_fname,1);
    indices = [];
    for jj = 1:N
        filename = [pathname char_fname(jj,:)];
        load(filename);
        indices = [indices; perform_computation(gl,char_fname(jj,:))];   
        tot_corr_trials = tot_corr_trials + gl.cum_correct_count(:,indices(end,:));
        tot_incorr_trials = tot_incorr_trials + gl.cum_incorrect_count(:,indices(end,:));
    end
%     keyboard;
    new_gl = gl;
    new_gl.cum_correct_count = tot_corr_trials;
    new_gl.cum_incorrect_count = tot_incorr_trials;
    new_gl.color_ID = gl.orig_color_ID;
    [~, alpha_all, beta_all] = perform_computation(new_gl,'Pooled Average');
    
else
    [indices, alpha_all, beta_all] = perform_computation(gl,'Pooled Average');
    tot_corr_trials = tot_corr_trials + gl.cum_correct_count(:,indices);
    tot_incorr_trials = tot_incorr_trials + gl.cum_incorrect_count(:,indices);
    new_gl = gl;
    new_gl.cum_correct_count = tot_corr_trials;
    new_gl.cum_incorrect_count = tot_incorr_trials;
    new_gl.color_ID = gl.orig_color_ID;
end

end

function [ind, alpha_all, beta_all] = perform_computation(gll,varargin)
figure_name = varargin{1};
new_idx = find_appro_index(gll);
prctcor = gll.cum_correct_count./(gll.cum_correct_count + gll.cum_incorrect_count);
prctcor = clearnans(prctcor);
new_x = logspace(log10(gll.stim.scale_lattice(1)),log10(gll.stim.scale_lattice(end)),101);
thresh = [];
model = [];
q = [];


% data and the psychometric function
ind = [];
for ii = 1:numel(gll.color_ID)
    ind = [ind find(new_idx==ii)];
    idx1 = ~isnan(prctcor(:,ind(end)));
    model_tmp = weibullFit_AD(gll.stim.scale_lattice(idx1)', [gll.cum_correct_count(idx1,ind(end)) gll.cum_incorrect_count(idx1,ind(end))], [],'single','fminsearch');
    model = [model; model_tmp];
    q = [q; 1 - 0.5.*exp(-((new_x./model(ii,1)).^model(ii,2)))];
    thresh = [thresh; model_tmp(1) model_tmp(2)];
end

mode = 0;
if mode == 1
    fig_h = figure('Visible','off');
    set(fig_h,'name',figure_name,'numbertitle','off')
    plot(gll.stim.scale_lattice, prctcor(:,ind(1)),'bo', 'MarkerFaceColor','b');hold on;
    plot(gll.stim.scale_lattice, prctcor(:,ind(2)),'mo', 'MarkerFaceColor','m');
    plot(gll.stim.scale_lattice, prctcor(:,ind(3)),'go', 'MarkerFaceColor','g');
    plot(gll.stim.scale_lattice, prctcor(:,ind(4)),'ko', 'MarkerFaceColor','k');
    plot(new_x', q(1,:),'b');
    plot(new_x', q(2,:),'m');
    plot(new_x', q(3,:),'g');
    plot(new_x', q(4,:),'k');
    xlabel('Contrast levels'), ylabel('Percentage correct choices');
    legend(gll.orig_color_ID{1},gll.orig_color_ID{2}, gll.orig_color_ID{3},gll.orig_color_ID{4} ,'location','SouthEast'); title('Psychometric Function');
    set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]);
    set(fig_h ,'Visible','on');
end
alpha_all = model(:,1);
beta_all = model(:,2);
tot_trials = gll.cum_correct_count + gll.cum_incorrect_count;

if strcmp(varargin{1},'Pooled Average')
    figure,subplot(221),plot(new_x', q(1,:),'b'); hold on; scatter(gll.stim.scale_lattice, prctcor(:,ind(1)),'bo','SizeData',tot_trials(:,ind(1))); 
    xlabel('Contrast'); ylabel('Per correct'); title(gll.orig_color_ID{1}); set(gca,'Xscale','log','Ylim',[0 1]);hold off;
    subplot(222),plot(new_x', q(2,:),'m'); hold on; scatter(gll.stim.scale_lattice, prctcor(:,ind(2)),'mo','SizeData',tot_trials(:,ind(2))); 
    xlabel('Contrast'); ylabel('Per correct'), title(gll.orig_color_ID{2}); set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    subplot(223),plot(new_x', q(3,:),'g'); hold on; scatter(gll.stim.scale_lattice, prctcor(:,ind(3)),'go', 'SizeData',tot_trials(:,ind(3))); 
    xlabel('Contrast'); ylabel('Per correct'), title(gll.orig_color_ID{3}); set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    subplot(224),plot(new_x', q(4,:),'k'); hold on; scatter(gll.stim.scale_lattice, prctcor(:,ind(4)),'ko', 'SizeData',tot_trials(:,ind(4))); 
    xlabel('Contrast'); ylabel('Per correct'), title(gll.orig_color_ID{4}); set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    
    % Just plotting the same result but now in a single figure
    figure,scatter(gll.stim.scale_lattice, prctcor(:,ind(1)),'bo','SizeData',tot_trials(:,ind(1))); hold on;
    scatter(gll.stim.scale_lattice, prctcor(:,ind(2)),'mo','SizeData',tot_trials(:,ind(2))); 
    scatter(gll.stim.scale_lattice, prctcor(:,ind(3)),'go', 'SizeData',tot_trials(:,ind(3))); 
    scatter(gll.stim.scale_lattice, prctcor(:,ind(4)),'ko', 'SizeData',tot_trials(:,ind(4))); 
    plot(new_x', q(1,:),'b');
    plot(new_x', q(2,:),'m');
    plot(new_x', q(3,:),'g');
    plot(new_x', q(4,:),'k');
    legend(gll.orig_color_ID{1},gll.orig_color_ID{2}, gll.orig_color_ID{3},gll.orig_color_ID{4} ,'location','SouthEast'); title('Psychometric Function');
    xlabel('Contrast'); ylabel('Percentage correct choices'), set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    
    % Plotting the the two poles of the same direction in the same plot,
    % basically one plot for orange-cyan and a second plot for lime-magenta
    % to see if there is any asymmetry
    % subplot 1 - orange-cyan
    % subplot 2 - lime-magenta
    % Also embedded is a statistical analysis (log-likelihood ratio test) to see how if the shift between the curves is statistically significant
%     keyboard;
    N = numel(gll.stim.scale_lattice);
    new_stim_lattice = repmat(gll.stim.scale_lattice',[2 1]);
    
    Indicator_lm = [ones(N,1); zeros(N,1)]; % 1- magenta, 0-lime
    lm_correct_count = [gll.cum_correct_count(:,ind(2)); gll.cum_correct_count(:,ind(4))]; 
    lm_incorrect_count = [gll.cum_incorrect_count(:,ind(2)); gll.cum_incorrect_count(:,ind(4))];
    lm_prctcor = lm_correct_count./(lm_correct_count + lm_incorrect_count);
    lm_trials = lm_correct_count + lm_incorrect_count;
    idxlm = ~isnan(lm_prctcor); 
    [modellm_com, fvallm_com,~,~] = weibullFit_AD(new_stim_lattice(idxlm), [lm_correct_count(idxlm) lm_incorrect_count(idxlm)], [], 'single', 'fminsearch');
    qlm_com = 1 - 0.5.*exp(-((new_x./modellm_com(1)).^modellm_com(2)));
    [modellm_ind, fvallm_ind,exitflaglm, hesslm] = weibullFit_AD(new_stim_lattice(idxlm), [lm_correct_count(idxlm) lm_incorrect_count(idxlm)], Indicator_lm(idxlm), 'multiple', 'fmincon',[modellm_com(1) modellm_com(2) modellm_com(1)]);
    qm = 1 - 0.5.*exp(-((new_x./modellm_ind(1)).^modellm_ind(2))); % fit for the lime data
    ql = 1 - 0.5.*exp(-((new_x./modellm_ind(3)).^modellm_ind(2))); % fit for the magenta data
    
    % The same thing for the orange cyan data
    Indicator_oc = [ones(N,1); zeros(N,1)]; % 1- orange, 0-cyan
    oc_correct_count = [gll.cum_correct_count(:,ind(1)); gll.cum_correct_count(:,ind(3))]; 
    oc_incorrect_count = [gll.cum_incorrect_count(:,ind(1)); gll.cum_incorrect_count(:,ind(3))];
    oc_prctcor = oc_correct_count./(oc_correct_count + oc_incorrect_count);
    oc_trials = oc_correct_count + oc_incorrect_count;
    idxoc = ~isnan(oc_prctcor); 
    [modeloc_com, fvaloc_com,~,~] = weibullFit_AD(new_stim_lattice(idxoc), [oc_correct_count(idxoc) oc_incorrect_count(idxoc)], [], 'single', 'fminsearch');
    qoc_com = 1 - 0.5.*exp(-((new_x./modeloc_com(1)).^modeloc_com(2)));
    [modeloc_ind, fvaloc_ind,exitflagoc, hessoc] = weibullFit_AD(new_stim_lattice(idxoc), [oc_correct_count(idxoc) oc_incorrect_count(idxoc)], Indicator_oc(idxoc), 'multiple', 'fmincon',[modeloc_com(1) modeloc_com(2) modeloc_com(1)]);
    qo = 1 - 0.5.*exp(-((new_x./modeloc_ind(1)).^modeloc_ind(2))); % fit for the lime data
    qc = 1 - 0.5.*exp(-((new_x./modeloc_ind(3)).^modeloc_ind(2))); % fit for the magenta data
    
    
    lm_LLratio = 2*(fvallm_com - fvallm_ind); % Degree of freedom
    oc_LLratio = 2*(fvaloc_com - fvaloc_ind); % Degree of freedom
    figure, subplot(221),plot(new_x', qlm_com,'b'); hold on; scatter(new_stim_lattice, lm_prctcor,'bo','SizeData',lm_trials); 
    title('Psychometric Function'); xlabel('Contrast'); ylabel('Percentage correct choices'), set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    subplot(222),plot(new_x', qoc_com,'b'); hold on;  scatter(new_stim_lattice, oc_prctcor,'bo','SizeData',oc_trials); 
    title('Psychometric Function'); xlabel('Contrast'); ylabel('Percentage correct choices'), set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    subplot(223),plot(new_x', qm,'m'); hold on; plot(new_x', ql,'k'); scatter(gll.stim.scale_lattice, lm_prctcor(1:N),'mo','SizeData',lm_trials(1:N));  
    scatter(gll.stim.scale_lattice, lm_prctcor(N+1:end),'ko','SizeData',lm_trials(N+1:end));  title('Psychometric Function'); xlabel('Contrast'); 
    ylabel('Percentage correct choices'); set(gca,'Xscale','log','Ylim',[0 1]); text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(lm_LLratio,1)))); 
    legend('magenta', 'lime' ,'location','SouthEast'); hold off;
    subplot(224),plot(new_x', qo,'b'); hold on; plot(new_x', qc,'g'); scatter(gll.stim.scale_lattice, oc_prctcor(1:N),'bo','SizeData',oc_trials(1:N));  
    scatter(gll.stim.scale_lattice, oc_prctcor(N+1:end),'go','SizeData',oc_trials(N+1:end));  title('Psychometric Function'); xlabel('Contrast'); 
    ylabel('Percentage correct choices'); set(gca,'Xscale','log','Ylim',[0 1]); text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(oc_LLratio,1)))); 
    legend('orange', 'cyan' ,'location','SouthEast'); hold off;
    
    
    % Trying to yoke lime and magenta together in order to determine a
    % combined threshold. Do the same for orange-cyan stimulus. I have also
    % embedded a log likelihood ratio test to see how significant is the
    % horizontal shift between the curves
    limemagenta_correct_trials = gll.cum_correct_count(:,ind(2)) + gll.cum_correct_count(:,ind(4));
    limemagenta_incorrect_trials = gll.cum_incorrect_count(:,ind(2)) + gll.cum_incorrect_count(:,ind(4));
    tottrials_limemagenta = limemagenta_correct_trials + limemagenta_incorrect_trials;
    prctcor_limemagenta = limemagenta_correct_trials./(limemagenta_correct_trials + limemagenta_incorrect_trials);
%     idx_lm = ~isnan(prctcor_limemagenta);
%     model_limemagenta = weibullFit_AD(gll.stim.scale_lattice(idx_lm)', [limemagenta_correct_trials(idx_lm) limemagenta_incorrect_trials(idx_lm)], [],'single','fminsearch');
%     modelfit_limemagenta = 1 - 0.5.*exp(-((new_x./model_limemagenta(1)).^model_limemagenta(2)));
    
    orangecyan_correct_trials = gll.cum_correct_count(:,ind(1)) + gll.cum_correct_count(:,ind(3));
    orangecyan_incorrect_trials = gll.cum_incorrect_count(:,ind(1)) + gll.cum_incorrect_count(:,ind(3));
    tottrials_orangecyan = orangecyan_correct_trials + orangecyan_incorrect_trials;
    prctcor_orangecyan = orangecyan_correct_trials./(orangecyan_correct_trials + orangecyan_incorrect_trials);
%     idx_oc = ~isnan(prctcor_orangecyan);
%     model_orangecyan = weibullFit_AD(gll.stim.scale_lattice(idx_oc)', [orangecyan_correct_trials(idx_oc) orangecyan_incorrect_trials(idx_oc)], [],'single','fminsearch');
%     modelfit_orangecyan = 1 - 0.5.*exp(-((new_x./model_orangecyan(1)).^model_orangecyan(2)));
    
 
    Indicator = [ones(N,1); zeros(N,1)]; % 1- orange-cyan, 2 - lime-magenta
    correct_count = [orangecyan_correct_trials; limemagenta_correct_trials];
    incorrect_count = [orangecyan_incorrect_trials; limemagenta_incorrect_trials];
    prctcor = correct_count./(correct_count + incorrect_count);
    trials = correct_count + incorrect_count;
    idx = ~isnan(prctcor); 
    [model_com, fval_com,~,~] = weibullFit_AD(new_stim_lattice(idx), [correct_count(idx) incorrect_count(idx)], [], 'single', 'fminsearch');
    q_com = 1 - 0.5.*exp(-((new_x./model_com(1)).^model_com(2)));
    [model_ind, fval_ind,exitflag, hess] = weibullFit_AD(new_stim_lattice(idx), [correct_count(idx) incorrect_count(idx)], Indicator(idx), 'multiple', 'fmincon');
    qoc = 1 - 0.5.*exp(-((new_x./model_ind(1)).^model_ind(2))); % fit for the orange-cyan data
    qlm = 1 - 0.5.*exp(-((new_x./model_ind(3)).^model_ind(2))); % fit for the lime-magenta data
    
    LLratio = 2*(fval_com - fval_ind);
    figure, subplot(121),plot(new_x', q_com,'b'); hold on; scatter(new_stim_lattice, prctcor,'bo','SizeData',trials); 
    title('Psychometric Function'); xlabel('Contrast'); ylabel('Percentage correct choices'), set(gca,'Xscale','log','Ylim',[0 1]); hold off; 
    subplot(122),plot(new_x', qoc,'b'); hold on; plot(new_x', qlm,'m'); scatter(gll.stim.scale_lattice, prctcor(1:N),'bo','SizeData',trials(1:N));  
    scatter(gll.stim.scale_lattice, prctcor(N+1:end),'mo','SizeData',trials(N+1:end));  title('Psychometric Function'); xlabel('Contrast'); 
    ylabel('Percentage correct choices'); set(gca,'Xscale','log','Ylim',[0 1]); text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(LLratio,1)))); 
    legend('orange-cyan', 'lime-magenta' ,'location','SouthEast'); hold off; 
      
    
    fitdata.modellm_ind = modellm_ind;
    fitdata.modeloc_ind = modeloc_ind;
    fitdata.model_ind = model_ind;
    modelErrs.lm = sqrt(abs(diag(inv(hesslm))));
    modelErrs.oc = sqrt(abs(diag(inv(hessoc))));
    modelErrs.model_ind = sqrt(abs(diag(inv(hessoc))));
    dataSLMTF_KM.fit_data = fitdata;
    dataSLMTF_KM.modelErrs = modelErrs;
    dataSLMTF_KM.gll = gll;
    dataSLMTF_KM.lm_LLratio = lm_LLratio;
    dataSLMTF_KM.oc_LLratio = oc_LLratio;
    dataSLMTF_KM.LLratio = LLratio;
    save('dataSLMTF_KM.mat','dataSLMTF_KM');
     
end

end