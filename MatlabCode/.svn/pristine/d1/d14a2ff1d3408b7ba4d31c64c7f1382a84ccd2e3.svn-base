function threshold = plot_fit(gll,fit_data)
% Abhishek - 05/2016
% This function is written with the purpose of just plotting the psychometric functions obtained from the fits within the function 'Analyze_data'.
% Only use the fits stored in the No_BACKUP directory


NCONTRASTS = 50; % 35 for ViewPixx data, 50 for Rig 2 Data
tot_corr_trials = zeros(NCONTRASTS,4);
tot_incorr_trials = zeros(NCONTRASTS,4);
if nargin < 2
    [fname2, pathname2] = uigetfile('*.mat', 'Select fit_data file');
    if isequal(fname2,0) || isequal(pathname2,0)
        return
    end
    char_fname = char(fname2);
    filename = [pathname2 char_fname];
    load(filename);
end
if nargin < 1
    [fname1, pathname1] = uigetfile('*.mat', 'Select gll structure');
    if isequal(fname1,0) || isequal(pathname1,0)
        return
    end
    char_fname = char(fname1);
    filename = [pathname1 char_fname];
    load(filename);
end

model1 = fit_data.model1; fval1 = fit_data.fval1;
model2 = fit_data.model2; fval2 = fit_data.fval2;
model3 = fit_data.model3; fval3 = fit_data.fval3;
model4 = fit_data.model4; fval4 = fit_data.fval4;

new_x = logspace(log10(gll.stim.scale_lattice(1)),log10(gll.stim.scale_lattice(end)),101);
Indicator = [ones(numel(gll.stim.scale_lattice),1); zeros(numel(gll.stim.scale_lattice),1)]; % 1 - L-M-S, 0 - L-M+S 
% Indicator = [ones(numel(gll.stim.scale_lattice),1); ones(numel(gll.stim.scale_lattice),1)];
N = numel(gll.stim.scale_lattice);

low_sf_correct_count = [gll.cum_correct_count(:,1); gll.cum_correct_count(:,2)];
low_sf_incorrect_count = [gll.cum_incorrect_count(:,1); gll.cum_incorrect_count(:,2)];
low_sf_prctcor = low_sf_correct_count./(low_sf_correct_count + low_sf_incorrect_count);
low_sf_trials = (low_sf_correct_count + low_sf_incorrect_count);

high_sf_correct_count = [gll.cum_correct_count(:,3); gll.cum_correct_count(:,4)];
high_sf_incorrect_count = [gll.cum_incorrect_count(:,3); gll.cum_incorrect_count(:,4)];
high_sf_prctcor = high_sf_correct_count./(high_sf_correct_count + high_sf_incorrect_count);
high_sf_trials = (high_sf_correct_count + high_sf_incorrect_count);


idx1 = ~isnan(low_sf_prctcor);
idx2 = ~isnan(high_sf_prctcor);
new_stim_lattice = repmat(gll.stim.scale_lattice',[2 1]);

% Obtaining the curve from the fitted parameters
q1 = 1 - 0.5.*exp(-((new_x./model1(1)).^model1(2))); % Low SF - Jointly fit the two directions
q2 = 1 - 0.5.*exp(-((new_x./model2(1)).^model2(2))); % High SF - Jointly fit the two directions
q31 = 1 - 0.5.*exp(-((new_x./model3(1)).^model3(2))); % Low SF - L-M-S
q32 = 1 - 0.5.*exp(-((new_x./model3(3)).^model3(2))); % Low SF - L-M+S
q41 = 1 - 0.5.*exp(-((new_x./model4(1)).^model4(2))); % High SF - L-M-S
q42 = 1 - 0.5.*exp(-((new_x./model4(3)).^model4(2))); % High SF - L-M+S

low_sf_LLratio = 2*(fval1 - fval3); % Degree of freedom
high_sf_LLratio = 2*(fval2 - fval4);
crit = chi2inv(0.95,1);
figure, axis square;
subplot(221), plot(new_x,q1); hold on; scatter(new_stim_lattice, low_sf_prctcor, 'bo', 'SizeData',low_sf_trials);
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices'); title ('Low sf'); hold off;
subplot(222), plot(new_x,q2); hold on; scatter(new_stim_lattice, high_sf_prctcor, 'bo', 'SizeData',high_sf_trials); 
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices'); title ('High sf'); hold off;
subplot(223), plot(new_x, q31,'b'); hold on; plot(new_x, q32,'m'); scatter(gll.stim.scale_lattice, low_sf_prctcor(1:N),'bo','SizeData',low_sf_trials(1:N)); scatter(gll.stim.scale_lattice, low_sf_prctcor(N+1:end) ,'mo','SizeData',low_sf_trials(N+1:end));
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices');  legend('L-M-S', 'L-M+S'); title ('Low sf');
text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(low_sf_LLratio,1)))); 
subplot(224), plot(new_x,q41,'b'); hold on; plot(new_x, q42,'m'); scatter(gll.stim.scale_lattice, high_sf_prctcor(1:N),'bo','SizeData',high_sf_trials(1:N)); scatter(gll.stim.scale_lattice, high_sf_prctcor(N+1:end) ,'mo','SizeData',high_sf_trials(N+1:end));
set(gca,'Xscale','log'); axis([gll.stim.scale_lattice(1) gll.stim.scale_lattice(end) 0 1]); xlabel('Contrast levels'), ylabel('Percentage correct choices'); legend('L-M-S', 'L-M+S'); title ('High sf'); 
text(0.3,0.1, strcat('p=',num2str(1-chi2cdf(high_sf_LLratio,1)))); hold off;
fprintf('low sf LL ratio is %1.2d with a p-value of %1.2d \n',low_sf_LLratio,1-chi2cdf(low_sf_LLratio,1));
fprintf('high sf LL ratio is %1.2d with a p-value of %1.2d \n',high_sf_LLratio,1-chi2cdf(high_sf_LLratio,1));
fprintf('95 percent confidence interval for chi-square dist = %1.2d \n',crit);

threshold.lowsf.soff = model3(1); % Low SF - L-M-S
threshold.lowsf.son = model3(3); % Low SF - L-M+S
threshold.highsf.soff = model4(1); % High SF - L-M-S
threshold.highsf.son = model4(3); % High SF - L-M+S


end

