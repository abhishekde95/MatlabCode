% Analyzing the contours of the firing rate surface from the Neurothresh files
% Author - Abhishek De, 10/18
% Not using this currently
close all; clearvars;
load filename_c.mat
load filename_l.mat
load THETA_all.mat
load RHO_all.mat 
load oog_idx_all.mat 
load not_oog_idx_all.mat
load linear_modelparams.mat
load quad_modelparams.mat
load S1LMS.mat
load S2LMS.mat
load WTS.mat 
load SpikeCounts.mat 
filename = [filename_c; filename_l];
linear_FRsurfacefitparams = [];
quad_FRsurfacefitparams = [];
linear_FRsurfacefitLL = [];
quad_FRsurfacefitLL = [];
diffmodelpredFR = [];
colormap('gray');
num_rows = 6;
count = 1;
plot_counter = 1;
for ii = 1:numel(filename)
    ind = ii;
    THETA = THETA_all{1,ind};
    THETA = THETA * pi/180; % converting to radians
    if any(THETA>(135*pi/180))
        allthetas = linspace(-pi,pi,100);
        newtheta = linspace(-pi,pi,101);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
        newtheta = linspace(-pi/4,3*pi/4,101);
    end
    RHO = RHO_all{1,ind};
    oog_idx = oog_idx_all{1,ind};
    not_oog_idx = not_oog_idx_all{1,ind};
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut);
    [x_orig, y_orig] = pol2cart(THETA,RHO);
    [~,eigval] = eig(cov([x_orig(not_oog_idx) y_orig(not_oog_idx)]));
    final_model1 = linear_modelparams(ind,:);
    final_model2 = quad_modelparams(ind,:);
                      
    % Need to fit the FR surfaces
    [linparams, fval1, FRsurface1,FRpredictedvalues1] = fitLinearmodeltospikecounts(cell2mat(WTS(1,ind)),cell2mat(SpikeCounts(1,ind)),final_model1);
    [quadparams, fval2, FRsurface2, FRpredictedvalues2] = fitQuadmodeltospikecounts(cell2mat(WTS(1,ind)),cell2mat(SpikeCounts(1,ind)),final_model2);
    linear_FRsurfacefitparams = [linear_FRsurfacefitparams; linparams];
    quad_FRsurfacefitparams = [quad_FRsurfacefitparams; quadparams];
    linear_FRsurfacefitLL = [linear_FRsurfacefitLL; fval1];
    quad_FRsurfacefitLL = [quad_FRsurfacefitLL; fval2];
    diffmodelpredFR = [diffmodelpredFR; mean((FRpredictedvalues1-FRpredictedvalues2).^2)];
    wts = cell2mat(WTS(1,ind));
    [X,Y] = meshgrid(linspace(-2,2,51));
    figure(plot_counter+1); 
    subplot(num_rows,3,3*count-2);surf(X,Y,FRsurface1); hold on; plot3(wts(:,1),wts(:,2),cell2mat(SpikeCounts(1,ind)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off;
    subplot(num_rows,3,3*count-1);surf(X,Y,FRsurface2); hold on; plot3(wts(:,1),wts(:,2),cell2mat(SpikeCounts(1,ind)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; 
    subplot(num_rows,3,3*count);fnplt(tpaps(wts',cell2mat(SpikeCounts(1,ind))')); hold on; plot3(wts(:,1),wts(:,2),cell2mat(SpikeCounts(1,ind)),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 1 0]); set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off;
    count  = count + 1;
    if count == (num_rows + 1)
        count = 1;
        plot_counter = plot_counter + 1;
    end
end

savevariables = 1;
if savevariables == 1
    save linear_FRsurfacefitparams linear_FRsurfacefitparams
    save quadFR_surfacefitparams quad_FRsurfacefitparams
    save linear_FRsurfacefitLL linear_FRsurfacefitLL
    save quad_FRsurfacefitLL quad_FRsurfacefitLL
    save diffmodelpredFR diffmodelpredFR
end