%% Function to load and plot MSE to compare OLS, LARS, and ASD

% Define names of datasets
expID = {'K021408004';...
    'K021508004';...
    'K031308004';...
    'K032108003'};

for iexpID = 1:length(expID)
    
    % Path to folder of results for each datapath
    ip2l = sprintf('D:/Users/Ryan/Desktop/%s_LARS_Figs/', expID{iexpID});
    
    % Load data from LARS and ASD
    LARSmse = load([ip2l 'LARS_MSE.mat']);
    ASDmse = load([ip2l 'ASD_MSE.mat']);
    
    % Ensure the OLS MSE is consistent between LARS and ASD
    if any(~ismembertol(LARSmse.STAmse, ASDmse.OLSmse))
        error('LARS and ASD did not return the same OLS mse');
    end
    
    % Plot all data
    figure; hold on;
    plot(LARSmse.frame_offset, LARSmse.STAmse, 'LineWidth',2);
    plot(LARSmse.frame_offset, LARSmse.CVmse, 'LineWidth', 2);
    plot(LARSmse.frame_offset, ASDmse.CVmse_aniso, 'LineWidth',2);
    hold off;
    legend('OLS STA','LARS', 'ASD');
    title([expID{iexpID} ' MSE'])
    ylabel('MSE'); xlabel('Frame Offset');
    saveas(gcf, [ip2l 'MSEcomp_LARS_ASD.png']);
end