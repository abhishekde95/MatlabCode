function Analyze_data(gl)
% close all;

% If u are reading this file, there is a chance that you might be confused
% which is a good analysis file one should refer for analysing SLMTF data (gl structure).
% You should use this file as it was written later
NCONTRASTS = 50; % 35 for ViewPixx data, 50 for Rig 2 Data and Propixx
tot_corr_trials = zeros(NCONTRASTS,4);
tot_incorr_trials = zeros(NCONTRASTS,4);
if nargin == 0
    [fname, pathname] = uigetfile('*.mat', 'Select a data file','MultiSelect','on');
    if isequal(fname,0) || isequal(pathname,0)
        return
    end
    char_fname = char(fname);
    N = size(char_fname,1);
    for jj = 1:N
        filename = [pathname char_fname(jj,:)];
        load(filename);
        indices = perform_computation(gl,char_fname(jj,:));   
        tot_corr_trials = tot_corr_trials + gl.cum_correct_count(:,indices);
        tot_incorr_trials = tot_incorr_trials + gl.cum_incorrect_count(:,indices);
    end
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

perform_loglikelihoodratio_test(new_gl, alpha_all, beta_all); % Compares a null and an alterate model based on the log likelihood ratio
end

function [sf_ind_tot, alpha_all, beta_all] = perform_computation(gll,varargin)

figure_name = varargin{1};
new_idx = find_appro_index(gll);
prctcor = gll.cum_correct_count./(gll.cum_correct_count + gll.cum_incorrect_count);
prctcor = clearnans(prctcor);
new_x = logspace(log10(gll.stim.scale_lattice(1)),log10(gll.stim.scale_lattice(end)),101);
N_sf = numel(gll.stim.sf);
thresh = [];
sf_ind_tot = [];
alpha_all = [];
beta_all = [];

for ii  = 1:N_sf
    sf_ind = [find(new_idx==N_sf*ii-1) find(new_idx==N_sf*ii)];
    
    % data and the psychometric function
    idx1 = ~isnan(prctcor(:,sf_ind(1)));
    model1 = weibullFit_AD(gll.stim.scale_lattice(idx1)', [gll.cum_correct_count(idx1,sf_ind(1)) gll.cum_incorrect_count(idx1,sf_ind(1))], [],'single','fminsearch');
    q1 = 1 - 0.5.*exp(-((new_x./model1(1)).^model1(2)));
    idx2 = ~isnan(prctcor(:,sf_ind(2)));
    model2 = weibullFit_AD(gll.stim.scale_lattice(idx2)', [gll.cum_correct_count(idx2,sf_ind(2)) gll.cum_incorrect_count(idx2,sf_ind(2))], [], 'single','fminsearch');
    q2 = 1 - 0.5.*exp(-((new_x./model2(1)).^model2(2)));
    alpha1 = model1(1); beta1 = model1(2);
    alpha2 = model2(1); beta2 = model2(2);
    thresh = [thresh; alpha1 alpha2];
    
    % graphical representation of the threshold in the 2 color directions
    theta = pi/4:pi/2:(2*pi)-pi/4;
    rho = zeros(size(theta));
    if gll.color_ID{sf_ind(1)} == strcat('(L-M)-S_sf',num2str(ii))
        rho(1) = alpha2;
        rho(3) = alpha2;
        rho(2) = alpha1;
        rho(4) = alpha1;
    elseif gll.color_ID{sf_ind(2)} == strcat('(L-M)+S_sf',num2str(ii))
        rho(1) = alpha1;
        rho(3) = alpha1;
        rho(2) = alpha2;
        rho(4) = alpha2;
    end
    
    theta = [theta, theta(1)];
    rho = [rho, rho(1)];
    
    radial_plot = 0;
    if radial_plot
        
        % plotting the figures
        fig_h = figure('Visible','off');
        set(fig_h,'name',figure_name,'numbertitle','off')
        
        subplot(121);plot(gll.stim.scale_lattice, prctcor(:,sf_ind(1)),'bo', 'MarkerFaceColor','b');hold on;
        plot(gll.stim.scale_lattice, prctcor(:,sf_ind(2)),'mo', 'MarkerFaceColor','m');
        plot(new_x', q1,'b');
        plot(new_x', q2,'m');
        xlabel('Contrast levels'), ylabel('Percentage correct choices');
        legend(gll.color_ID{sf_ind(1)}, gll.color_ID{sf_ind(2)},'location','SouthEast'); title('Psychometric Function');
        set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]);
        
        subplot(122);
        polar(theta,rho,'k-'); hold off;ylabel('S'),xlabel('L-M');
        title('Threshold plot');hold off;
        
        % Useful for writing texts in the figure
        uicontrol(fig_h,'Position',[320 350 150 20],'Style','text',...
            'String','Spatial_frequency - ','Units','pixels');
        uicontrol(fig_h,'Position',[440 350 20 20],'Style','text',...
            'String',num2str(gll.stim.sf(ii)),'Units','pixels');
        uicontrol(fig_h,'Position',[330 50 80 20],'Style','text',...
            'String','Threshold 1 - ','Units','pixels');
        uicontrol(fig_h,'Position',[400 50 50 20],'Style','text',...
            'String',num2str(thresh(ii,1)),'Units','pixels');
        uicontrol(fig_h,'Position',[330 30 80 20],'Style','text',...
            'String','Threshold 2 - ','Units','pixels');
        uicontrol(fig_h,'Position',[400 30 50 20],'Style','text',...
            'String',num2str(thresh(ii,2)),'Units','pixels');
        set(fig_h ,'Visible','on');
    end
    sf_ind_tot = [sf_ind_tot sf_ind];
    alpha_all = [alpha_all; alpha1 alpha2];
    beta_all = [beta_all; beta1 beta2];
end
end


function perform_loglikelihoodratio_test(gll,alpha_guess, beta_guess)
new_x = logspace(log10(gll.stim.scale_lattice(1)),log10(gll.stim.scale_lattice(end)),101);
Indicator = [ones(numel(gll.stim.scale_lattice),1); zeros(numel(gll.stim.scale_lattice),1)]; % 1 - L-M-S, 0 - L-M+S 
% Indicator = [ones(numel(gll.stim.scale_lattice),1); ones(numel(gll.stim.scale_lattice),1)];
N = numel(gll.stim.scale_lattice);

low_sf_correct_count = [gll.cum_correct_count(:,1); gll.cum_correct_count(:,2)];
low_sf_incorrect_count = [gll.cum_incorrect_count(:,1); gll.cum_incorrect_count(:,2)];
low_sf_prctcor = low_sf_correct_count./(low_sf_correct_count + low_sf_incorrect_count);
low_sf_trials = low_sf_correct_count + low_sf_incorrect_count;

high_sf_correct_count = [gll.cum_correct_count(:,3); gll.cum_correct_count(:,4)];
high_sf_incorrect_count = [gll.cum_incorrect_count(:,3); gll.cum_incorrect_count(:,4)];
high_sf_prctcor = high_sf_correct_count./(high_sf_correct_count + high_sf_incorrect_count);
high_sf_trials = high_sf_correct_count + high_sf_incorrect_count;

idx1 = ~isnan(low_sf_prctcor);
idx2 = ~isnan(high_sf_prctcor);
new_stim_lattice = repmat(gll.stim.scale_lattice',[2 1]);

% Fitting the data: no distinction between the two datasets
[model1, fval1,~,~] = weibullFit_AD(new_stim_lattice(idx1), [low_sf_correct_count(idx1) low_sf_incorrect_count(idx1)], [], 'single', 'fminsearch', [max(alpha_guess(1,:)), min(beta_guess(1,:))]);
q1 = 1 - 0.5.*exp(-((new_x./model1(1)).^model1(2)));
[model2, fval2,~, ~] = weibullFit_AD(new_stim_lattice(idx2), [high_sf_correct_count(idx2) high_sf_incorrect_count(idx2)], [],'single', 'fminsearch', [alpha_guess(2,2), beta_guess(2,2)]);
q2 = 1 - 0.5.*exp(-((new_x./model2(1)).^model2(2)));

% Fitting the data individually for each direction using an indicator function, might be some kind of error in weibull_AD
% Fitting the low sf data
[model3, fval3,exitflag3, hess3] = weibullFit_AD(new_stim_lattice(idx1), [low_sf_correct_count(idx1) low_sf_incorrect_count(idx1)], Indicator(idx1), 'multiple', 'fmincon');%, [model2(1), model2(2), model2(1)]);
q31 = 1 - 0.5.*exp(-((new_x./model3(1)).^model3(2)));
q32 = 1 - 0.5.*exp(-((new_x./model3(3)).^model3(2)));

% Fitting the high sf data
[model4, fval4, exitflag4, hess4] = weibullFit_AD(new_stim_lattice(idx2), [high_sf_correct_count(idx2) high_sf_incorrect_count(idx2)], Indicator(idx2), 'multiple', 'fmincon');%, [model2(1), model2(2), model2(1)]);
q41 = 1 - 0.5.*exp(-((new_x./model4(1)).^model4(2)));
q42 = 1 - 0.5.*exp(-((new_x./model4(3)).^model4(2)));
low_sf_LLratio = 2*(fval1 - fval3); % Degree of freedom
high_sf_LLratio = 2*(fval2 - fval4);


figure, subplot(221), plot(new_x,q1); hold on; scatter(new_stim_lattice, low_sf_prctcor, 'bo','SizeData',low_sf_trials);
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices'); title ('Low sf'); hold off;
subplot(222), plot(new_x,q2); hold on; scatter(new_stim_lattice, high_sf_prctcor, 'bo','SizeData',high_sf_trials); 
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices'); title ('High sf'); hold off;

subplot(223), plot(new_x, q31,'b'); hold on; plot(new_x, q32,'m'); scatter(gll.stim.scale_lattice, low_sf_prctcor(1:N),'bo','SizeData',low_sf_trials(1:N)); scatter(gll.stim.scale_lattice, low_sf_prctcor(N+1:end) ,'mo','SizeData',low_sf_trials(N+1:end));
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices'); title ('Low sf'); text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(low_sf_LLratio,1)))); hold off;
subplot(224), plot(new_x,q41,'b'); hold on; plot(new_x, q42,'m'); scatter(gll.stim.scale_lattice, high_sf_prctcor(1:N),'bo','SizeData',high_sf_trials(1:N)); scatter(gll.stim.scale_lattice, high_sf_prctcor(N+1:end) ,'mo','SizeData',high_sf_trials(N+1:end));
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices'); title ('High sf'); text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(high_sf_LLratio,1)))); hold off;


crit = chi2inv(0.95,1);

fprintf('low sf LL ratio is %1.2d \n',low_sf_LLratio);
fprintf('high sf LL ratio is %1.2d \n',high_sf_LLratio);
fprintf('95 percent confidence interval for chi-square dist = %1.2d \n',crit);

% saving the parameters
fit_data.model1 = model1; fit_data.fval1 = fval1;
fit_data.model2 = model2; fit_data.fval2 = fval2;
fit_data.model3 = model3; fit_data.fval3 = fval3;
fit_data.model4 = model4; fit_data.fval4 = fval4;
modelErrs.low_sf = sqrt(abs(diag(inv(hess3))));
modelErrs.high_sf = sqrt(abs(diag(inv(hess4))));
dataSLMTF.fit_data = fit_data;
dataSLMTF.modelErrs = modelErrs;
dataSLMTF.gll = gll;
dataSLMTF.low_sf_LLratio = low_sf_LLratio;
dataSLMTF.high_sf_LLratio = high_sf_LLratio;
save('dataSLMTF.mat','dataSLMTF');
end


