function Analyze_data_KMwsf(gl)
 clearvars;
% This is an extension of the SLMTF_KM where we are now extracting the
% detection thresholds for 0.5 cycles/deg orange-cyan and lime-magenta
% modulations
NCONTRASTS = 50; % 50 for Rig 1 Data and Propixx
tot_corr_trials = zeros(NCONTRASTS,2);
tot_incorr_trials = zeros(NCONTRASTS,2);
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
    [indices, alpha_all, beta_all] = perform_computation(gl,'');
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
    plot(new_x', q(1,:),'b');
    plot(new_x', q(2,:),'m');
    xlabel('Contrast levels'), ylabel('Percentage correct choices');
    legend(gll.orig_color_ID{1},gll.orig_color_ID{2},'location','SouthEast'); title('Psychometric Function');
    set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]);
    set(fig_h ,'Visible','on');
end
alpha_all = model(:,1);
beta_all = model(:,2);
tot_trials = gll.cum_correct_count + gll.cum_incorrect_count;

if strcmp(varargin{1},'Pooled Average')
    figure,subplot(121),plot(new_x', q(1,:),'b'); hold on; scatter(gll.stim.scale_lattice, prctcor(:,ind(1)),'bo','SizeData',tot_trials(:,ind(1))); 
    xlabel('Contrast'); ylabel('Per correct'); title(gll.orig_color_ID{1}); set(gca,'Xscale','log','Ylim',[0 1]);hold off;
    subplot(122),plot(new_x', q(2,:),'m'); hold on; scatter(gll.stim.scale_lattice, prctcor(:,ind(2)),'mo','SizeData',tot_trials(:,ind(2))); 
    xlabel('Contrast'); ylabel('Per correct'), title(gll.orig_color_ID{2}); set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    
    % Just plotting the same result but now in a single figure
    figure,scatter(gll.stim.scale_lattice, prctcor(:,ind(1)),'bo','SizeData',tot_trials(:,ind(1))); hold on;
    scatter(gll.stim.scale_lattice, prctcor(:,ind(2)),'mo','SizeData',tot_trials(:,ind(2))); 
    plot(new_x', q(1,:),'b');
    plot(new_x', q(2,:),'m');
    legend(gll.orig_color_ID{1},gll.orig_color_ID{2}, 'location','SouthEast'); title('Psychometric Function');
    xlabel('Contrast'); ylabel('Percentage correct choices'), set(gca,'Xscale','log','Ylim',[0 1]); hold off;
    
    % Plotting the the two poles of the same direction in the same plot,
    % basically one plot for orange-cyan and a second plot for lime-magenta
    % to see if there is any asymmetry
    % subplot 1 - orange-cyan
    % subplot 2 - lime-magenta
    % Also embedded is a statistical analysis (log-likelihood ratio test) to see how if the shift between the curves is statistically significant
    N = numel(gll.stim.scale_lattice);
    new_stim_lattice = repmat(gll.stim.scale_lattice',[2 1]);
    
    Indicator = [ones(N,1); zeros(N,1)]; % 1- orange-cyan, 0- lime-magenta
    correct_count_com = [gll.cum_correct_count(:,ind(1)); gll.cum_correct_count(:,ind(2))]; 
    incorrect_count_com = [gll.cum_incorrect_count(:,ind(1)); gll.cum_incorrect_count(:,ind(2))];
    prctcor_com = correct_count_com./(correct_count_com + incorrect_count_com);
    tottrials = correct_count_com + incorrect_count_com;
    idx = ~isnan(prctcor_com); 
    [model_com, fval_com,~,~] = weibullFit_AD(new_stim_lattice(idx), [correct_count_com(idx) incorrect_count_com(idx)], [], 'single', 'fminsearch');
    q_com = 1 - 0.5.*exp(-((new_x./model_com(1)).^model_com(2)));
    [model_ind, fval_ind,exitflag, hess] = weibullFit_AD(new_stim_lattice(idx), [correct_count_com(idx) incorrect_count_com(idx)], Indicator(idx), 'multiple', 'fminsearch');
    qoc = 1 - 0.5.*exp(-((new_x./model_ind(1)).^model_ind(2))); % fit for the orange-cyan data
    qlm = 1 - 0.5.*exp(-((new_x./model_ind(3)).^model_ind(2))); % fit for the lime-magenta data
    
    LLratio = 2*(fval_com - fval_ind); % Degree of freedom
    figure, subplot(121),plot(new_x', q_com,'b'); hold on; scatter(new_stim_lattice, prctcor_com,'bo','SizeData',tottrials); 
    title('Psychometric Function'); xlabel('Contrast'); ylabel('Percentage correct choices'), set(gca,'Xscale','log','Ylim',[0 1]); hold off;   
    subplot(122),plot(new_x', qoc,'b'); hold on; plot(new_x', qlm,'m'); scatter(gll.stim.scale_lattice, prctcor_com(1:N),'bo','SizeData',tottrials(1:N));  
    scatter(gll.stim.scale_lattice, prctcor_com(N+1:end),'mo','SizeData',tottrials(N+1:end));  title('Psychometric Function'); xlabel('Contrast'); 
    ylabel('Percentage correct choices'); set(gca,'Xscale','log','Ylim',[0 1]); text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(LLratio,1)))); 
    legend('orange-cyan', 'lime-magenta' ,'location','SouthEast'); hold off;
    
    fitdata.model_ind = model_ind;
    modelErrs = sqrt(abs(diag(inv(hess))));
    dataSLMTF_KMwsf.fit_data = fitdata;
    dataSLMTF_KMwsf.modelErrs = modelErrs;
    dataSLMTF_KMwsf.gll = gll;
    dataSLMTF_KMwsf.LLratio = LLratio;
    save('dataSLMTF_KMwsf.mat','dataSLMTF_KMwsf');
     
end

end
