% This script will allow some model predictions on FA laser vs FA control
% trials based on some ROC analysis
% Author - Abhishek De, 6/18
close all;
clearvars;
mus = repmat([-0.5 0 0.5],[3 1]);
sigmas = repmat([0.5; 1; 2],[1 3]);
plot_counter = 1;
xvals = linspace(-4,4,101);
for ii = 1:3
    for jj = 1:3
        figure(plot_counter); subplot(3,3,3*(ii-1)+jj); plot(xvals,normpdf(xvals,mus(ii,jj),sigmas(ii,jj)),'g','Linewidth',2); hold on; 
        plot(xvals,normpdf(xvals,0,1),'k','Linewidth',2); axis square; title(strcat('mu=',num2str(mus(ii,jj)),',sig=',num2str(sigmas(ii,jj))));hold off;
        
        figure(plot_counter+1); subplot(3,3,3*(ii-1)+jj); plot(1-normcdf(xvals,mus(ii,jj),sigmas(ii,jj)),1-normcdf(xvals,0,1),'b','Linewidth',2); hold on;
        line([0 1],[0 1],'Color','k','LineStyle','--'); axis square; title(strcat('mu=',num2str(mus(ii,jj)),',sig=',num2str(sigmas(ii,jj))));hold off; 
        
    end
end
plot_counter = plot_counter + 2;

% Now I want to see if any tranformation to the the two distributions affect the curve or not
for ii = 1:3
    for jj = 1:3
        figure(plot_counter); subplot(2,3,1);plot(xvals,normpdf(xvals,mus(ii,jj),sigmas(ii,jj)),'g','Linewidth',2); hold on;
        plot(xvals,normpdf(xvals,0,1),'k','Linewidth',2); axis square; title(strcat('mu=',num2str(mus(ii,jj)),',sig=',num2str(sigmas(ii,jj))));hold off;
        subplot(2,3,2);plot(xvals,normpdf(xvals/3,mus(ii,jj),sigmas(ii,jj)),'g','Linewidth',2); hold on; plot(xvals,normpdf(xvals/3,0,1),'k','Linewidth',2); title('stretched'); axis square; hold off;
        subplot(2,3,3);plot(xvals,normpdf(xvals*3,mus(ii,jj),sigmas(ii,jj)),'g','Linewidth',2); hold on; plot(xvals,normpdf(xvals*3,0,1),'k','Linewidth',2); title('contract'); axis square; hold off;
        subplot(2,3,4); plot(1-normcdf(xvals,mus(ii,jj),sigmas(ii,jj)),1-normcdf(xvals,0,1),'b','Linewidth',2); hold on;
        line([0 1],[0 1],'Color','k','LineStyle','--'); axis square; title(strcat('mu=',num2str(mus(ii,jj)),',sig=',num2str(sigmas(ii,jj))));hold off;
        subplot(2,3,5); plot(1-normcdf(xvals/3,mus(ii,jj),sigmas(ii,jj)),1-normcdf(xvals/3,0,1),'b','Linewidth',2); hold on; line([0 1],[0 1],'Color','k','LineStyle','--'); title('stretch'); axis square; hold off;
        subplot(2,3,6); plot(1-normcdf(xvals*3,mus(ii,jj),sigmas(ii,jj)),1-normcdf(xvals*3,0,1),'b','Linewidth',2); hold on; line([0 1],[0 1],'Color','k','LineStyle','--'); title('contract'); axis square; hold off;
        plot_counter = plot_counter + 1;
    end
end
% The curves are invariant to changes in the distribution - linear
% stretching and contraction