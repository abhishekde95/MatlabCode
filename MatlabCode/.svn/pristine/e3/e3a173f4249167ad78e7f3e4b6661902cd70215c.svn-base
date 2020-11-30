npts = 20;
support_x = 1:npts;
support_y = 1:npts;
[xx, yy] = meshgrid(support_x, support_y);

support = [xx(:) yy(:)];

% want to make a somewhat spiky function to cause numerical problem
muvec1 = [0.4*mean(support_x) 0.3*mean(support_y)];
muvec2 = [0.9*mean(support_x) 1*mean(support_y)];
muvec3 = [1.5*mean(support_x) 1.6*mean(support_y)];

invC1 = diag([mean(support_x)/40; mean(support_y)/40]);
invC2 = diag([mean(support_x)/40; mean(support_y)/60]);
invC3 = diag([mean(support_x)/50; mean(support_y)/30]);

const = 0.1;
intensity = @(t) const + 15*diag(exp(-0.5*(bsxfun(@minus, t, muvec1))*invC1*(bsxfun(@minus, t, muvec1))')) + 30*diag(exp(-0.5*(bsxfun(@minus, t, muvec2))*invC2*(bsxfun(@minus, t, muvec2))')) + 20*diag(exp(-0.5*(bsxfun(@minus, t, muvec3))*invC3*(bsxfun(@minus, t, muvec3))'));
lambda_true =intensity(support);

subplot(241);
imagesc(support_x, support_y, reshape(lambda_true, npts, []));axis image; axis xy;
title('true lambda');

%%

clear;
clc;

maxrep = 100;
HowmanyData = 100;

load mse_mi;
mse_mi_logexp = mean(mse_mi,2);

load ../NewAlgorithm_NoKinv_ExpG/mse_mi;
mse_mi_exp = mean(mse_mi);

load mse_br;
mse_br_logexp = mean(mse_br,2);

load ../NewAlgorithm_NoKinv_ExpG/mse_br;
mse_br_exp = mean(mse_br);

subplot(242);
plot(mse_mi_logexp, 'b.');title('mse');
hold on;
plot(mse_br_logexp, 'r.');
set(gca, 'ylim', [20 110]);

subplot(246);
plot(mse_mi_exp, 'b.');
hold on;
plot(mse_br_exp, 'r.');
set(gca, 'ylim', [20 110]);


%%
load totVar_mi;
totVar_mi_logexp = mean(totVar_mi,2);

load ../NewAlgorithm_NoKinv_ExpG/totVar_mi;
totVar_mi_exp = mean(totVar_mi);

load totVar_br;
totVar_br_logexp = mean(totVar_br,2);

load ../NewAlgorithm_NoKinv_ExpG/totVar_br;
totVar_br_exp = mean(totVar_br);

subplot(243);
semilogy(totVar_mi_logexp, 'b.');title('tot var');
hold on;
semilogy(totVar_br_logexp, 'r.'); grid on;
set(gca, 'ylim',[10^2 1.2*10^4]);

subplot(247);
semilogy(totVar_mi_exp, 'b.');
hold on;
semilogy(totVar_br_exp, 'r.'); grid on;
set(gca, 'ylim',[10^2 1.2*10^4]);

%%
load postent_mi;
postent_mi_logexp = mean(postent_mi,2);

load ../NewAlgorithm_NoKinv_ExpG/postent_mi;
postent_mi_exp = mean(postent_mi);

load postent_br;
postent_br_logexp = mean(postent_br,2);

load ../NewAlgorithm_NoKinv_ExpG/postent_br;
postent_br_exp = mean(postent_br);

subplot(244);
plot(postent_mi_logexp, 'b.'); 
hold on;
plot(postent_br_logexp, 'r.');

subplot(248);
plot(postent_mi_exp, 'b.');
hold on;
plot(postent_br_exp, 'r.');
set(gca,'ylim',[-6350 -6210]);

%%

figure(1);
subplot(222);
plot(1:HowmanyData, mean(mse_mi,2), '-b.', 1:HowmanyData, mean(mse_br,2), '-r.');
% set(gca, 'ylim', [40 180]);
title('mse'); legend('mutual info', 'bayes risk (BR)');
xlabel('# samps');

% load totVar_us;
load totVar_mi;
load totVar_br;

subplot(223);
semilogy(1:HowmanyData, mean(totVar_mi,2), '-b.', 1:HowmanyData, mean(totVar_br,2), '-r.');
title('total variance');
% set(gca, 'ylim', [650 10^5]);
% legend('uncertainty sampling', 'mutual info', 'bayes risk (BR)');
% xlabel('# samps');


% load postent_us;
load postent_mi;
load postent_br;

subplot(224);
plot( 1:HowmanyData, mean(postent_mi,2), '-b.', 1:HowmanyData, mean(postent_br,2), '-r.');
title('posterior entropy');
% set(gca, 'ylim', [35 255]);
% title('mse'); legend('uncertainty sampling', 'mutual info', 'bayes risk (BR)');
xlabel('# samps');