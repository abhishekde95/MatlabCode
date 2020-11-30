% Writing a script to test sensitivity vs criterion change
% Author - Abhishek De, 6/19

close all; clearvars;
plot_counter = 1;
mu1 = 0;
mu2 = 1;
sigma = 1;

x = linspace(-4,4,101);
dprime1 = []; dprime2  = []; dprime3 = []; dprime4 = [];
C1 = []; C2 = []; C3 = []; C4 = [];
for ii = -3:0.2:3.0

    H1 = 1-normcdf(ii,mu2,sigma);
    FA1 = 1-normcdf(ii,mu1,sigma);
    dprime1 = [dprime1; norminv(H1)-norminv(FA1)];
    C1 = [C1; 0.5*(norminv(H1)+norminv(FA1))];
    
    H2 = 1-normcdf(0.5,ii,sigma);
    FA2 = 1-normcdf(0.5,mu1,sigma);
    dprime2 = [dprime2; norminv(H2)-norminv(FA2)];
    C2 = [C2; 0.5*(norminv(H2)+norminv(FA2))];
    
    H3 = 1-normcdf(0,ii,sigma);
    FA3 = 1-normcdf(0,mu1,sigma);
    dprime3 = [dprime3; norminv(H3)-norminv(FA3)];
    C3 = [C3; 0.5*(norminv(H3)+norminv(FA3))];
    
    H4 = 1-normcdf(1,ii,sigma);
    FA4 = 1-normcdf(1,mu1,sigma);
    dprime4 = [dprime4; norminv(H4)-norminv(FA4)];
    C4 = [C4; 0.5*(norminv(H4)+norminv(FA4))];
    
end

figure(plot_counter);set(gcf,'Name','Changing criterion & sensitivity');
subplot(121); plot(dprime1,C1,'Color',[0 0 0],'MarkerEdgeColor',[1 1 1]); xlabel('d prime1'); ylabel('C1'); title('Changing criterion'); set(gca,'Tickdir','out','Xlim',[-3 3],'Ylim',[-3 3]); axis square;
subplot(122); plot(dprime2,C2,'Color',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(dprime3,C3,'Color',[0 1 0],'MarkerEdgeColor',[1 1 1]); plot(dprime4,C4,'Color',[0 0 1],'MarkerEdgeColor',[1 1 1]);
xlabel('d prime'); ylabel('Crit'); legend('crit = 0.5','crit = 0','crit = 1.0'); title('Changing d prime'); set(gca,'Tickdir','out','Xlim',[-3 3],'Ylim',[-3 3]); axis square;
plot_counter = plot_counter + 1;
