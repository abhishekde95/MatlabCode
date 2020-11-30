%%
% Trying to simulate (and then fit) Abhishek, neurothresh-type data 
% First generating fake data
weights = [1 10];
GAMUTEDGE = 20; 
thetas = linspace(-pi/4, 3*pi/4,20);
error_var = 0.2;
dotprods = weights*[cos(thetas); sin(thetas)];
staircase_terminations_mn = 1./dotprods;
% Sanity check: Once I scale these unit vectors by
% "staircase_terminations_mn" they should all have the same dot product
% onto "weights". staircase_terminations_mn is the expected termination 
% distance for each of the adaptive searches.
% weights*(repmat(staircase_terminations_mn,2,1).*[cos(thetas); sin(thetas)]);
% Some searches go out of gamut
L =staircase_terminations_mn > GAMUTEDGE | staircase_terminations_mn < 0;
LOOG = false(length(thetas),1);
LOOG(L) = true;
% adding some (log-normal) error
mu = log(staircase_terminations_mn(~LOOG).^2./sqrt(error_var+staircase_terminations_mn(~LOOG).^2));
sigma = sqrt(log(error_var./staircase_terminations_mn(~LOOG).^2 + 1));
staircase_terminations = zeros(length(thetas),1);
staircase_terminations(LOOG) = GAMUTEDGE;
staircase_terminations(~LOOG) = lognrnd(mu,sigma);
% At this point "staircase_terminations" is the distance of each of the
% searches (each in a direction theta)
% Here's some plotting stuff that I commented out that plots the fake data
% in cartisian and polar coordinates.
figure; subplot(2,1,1); hold on; axis square;
[x,y] = pol2cart(thetas', staircase_terminations);
plot(x(~LOOG),y(~LOOG),'ko');
plot(x(LOOG),y(LOOG),'ro');
% 
% subplot(2,1,2); hold on; axis square;
% plot(thetas(~LOOG),log10(staircase_terminations(~LOOG)),'ko');
% plot(thetas(LOOG),log10(staircase_terminations(LOOG)),'ro');
% Calculating a goodness of fit for a set of model parameters (grid search)
v = linspace(-20,20,100);
data = zeros(length(v),length(v));
for i = 1:length(v)
    for j = 1:length(v)
        modelparams = [v(i) v(j)];
        pred_staircase_terminations = 1./(modelparams*[cos(thetas); sin(thetas)])';
        L = pred_staircase_terminations > GAMUTEDGE | pred_staircase_terminations < 0;
        pred_staircase_terminations(L) = GAMUTEDGE;
        % Ignoring OOG points in the calculation of error for now.
        if sum(L)./length(L) > 0.5 % if > 0.5 of the seaches go out of gamut, this is probably not a good fit (e.g. weights are too small)
            err = nan;
        else
            err = mean((log(pred_staircase_terminations(~L))-log(staircase_terminations(~L))).^2);
        end
        data(i,j) = err;
    end
end
figure; 
surf(data);
[tmp_i, tmp_j] = ind2sub([length(v) length(v)],find(data==min(min(data))));
bestparams = [v(tmp_i) v(tmp_j)];
figure; subplot(2,1,1); hold on; axis square;
[x,y] = pol2cart(thetas', staircase_terminations);
plot(x(~LOOG),y(~LOOG),'ko');
plot(x(LOOG),y(LOOG),'ro');
allthetas = linspace(-pi,pi,100);
preds = 1./(bestparams*[cos(allthetas); sin(allthetas)]);
LOOGtmp= preds>GAMUTEDGE|preds<0;
[x,y] = pol2cart(allthetas(~LOOGtmp)', preds(~LOOGtmp)');
plot(x,y,'b-');
plot(0,0,'m*');
subplot(2,1,2); hold on; axis square;
plot(thetas(~LOOG),log10(staircase_terminations(~LOOG)),'ko');
plot(thetas(LOOG),log10(staircase_terminations(LOOG)),'ro');
plot(allthetas(~LOOGtmp)', log10(preds(~LOOGtmp)'),'b-');
weights
bestparams
% Hmmm... seems to be working OK, but bestparams is definitely biased.