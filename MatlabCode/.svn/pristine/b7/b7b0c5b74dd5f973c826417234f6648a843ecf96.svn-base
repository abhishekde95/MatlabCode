% Script for testing bimodality using a bootstrap method: based on Efron
% and Tibshirani's method. Also implementing the cross-validation procedure
% on bin size based on Rich Pang's suggestion
% Author - Abhishek De, 11/19

close all; clearvars;
load Output_List_waveform.mat
load timediff.mat
plot_counter = 1;
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
timediff(~Z_cellsofinterest) = [];

% The first part is about determining the optimal bin-width for the kernel density function
K = 10;
c = cvpartition(timediff,'KFold',K);
h = 5:5:100;
p = zeros(numel(h),1);
cost = [];
for kk = 1:numel(h)
    lik = [];
    for ii = 1:K
        vec1 = 10:h(kk):700;
        vec = vec1(1:end-1)+vec1(2:end);
        L = numel(vec);
        x_train = timediff(c.training(ii));
        x_test = timediff(c.test(ii));
        N = numel(x_train);
        d_train = zeros(L,1);
        for jj = 1:L
            d_train(jj) = (1/(N*h(kk)))*sum(exp((-((vec(jj)-x_train)/h(kk)).^2/2)*sqrt(2)*pi));
        end
        d_train(d_train==0) = eps;
        N = histcounts(x_test,vec1);    
        lik = [lik; mean(-log(d_train).*N')]; % negative log-likelihood 
        
    end
    cost = [cost; mean(lik)];
end
[~,I] = max(cost);
h_opt = h(I);
figure(1); subplot(131); plot(h,cost,'k','Linewidth',2); axis square; set(gca,'Tickdir','out'); xlabel('kernel width'); ylabel('neg log lik');
subplot(132); histogram(timediff,10:h_opt:700,'FaceColor','k','EdgeColor',[1 1 1]); axis square; set(gca,'Tickdir','out','Xlim',[0 700]); ylabel('count'); xlabel('waveform duration');

% Now testing whether the probability 
nboot = 500;
[p,occurences] = testbimodality(h_opt,timediff,10,700,nboot);

pval = logspace(-20,-1,51);
y = nboot*binopdf(0,nboot,pval);
p_opt = pval(find(y<occurences(2),1));
subplot(133); plot(pval,y,'k','Linewidth',2); hold on; plot(p_opt,0,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','XScale','log'); xlabel('p value'); ylabel('Counts of bimodality');  title('Binomial dist');
plot_counter = plot_counter + 1;


