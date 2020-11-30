%used in LMTFPaperFigures - not used to generate a figure that was used in
%the final paper, but probably still important to have around.
function gaussianProcessRegressionFitLMTF(data, zoom, bkgndrgb, M)
%taken from GBS Section 17 and rawDataPlot
zoom_data = [0,0,0,0,0,0,0];
covfunc = {'covProd',  {{'covMask',{[1 0],'covSEiso'}}, {'covMask',{[0 1],'covPeriodPi'}}}};
meanfunc = @meanConst;
likfunc = @likGauss;

hyp.cov = [3;0;1;-1];
hyp.mean = .05;
hyp.lik = -1;
if zoom
    for i = 1:length(data)
        if data(i,3) <= 10.000
            zoom_data = [zoom_data; data(i,:)];
        else
            continue
        end
    end
    zoom_data(1,:) = [];
    data = zoom_data;
end

Loog = logical(data(:,4));

% Plotting the data in LMTF space
figure; axes; hold on;
plot3(data(~Loog,1),data(~Loog,2),data(~Loog,3),'k.');
plot3(-data(~Loog,1),-data(~Loog,2),data(~Loog,3),'k.');
%plot3(data(Loog,1),data(Loog,2),data(Loog,3),'k.');
%plot3(-data(Loog,1),-data(Loog,2),data(Loog,3),'k.');

% Plotting the final regression fit in LMTF space
[xx,yy,zz] = meshgrid(linspace(-max(abs(data(:,1))),max(abs(data(:,1))),size(data,1)),...
    linspace(-max(abs(data(:,2))),max(abs(data(:,2))),size(data,1)),...
    linspace(min(data(:,3)),max(data(:,3)),size(data,1)));
[th_star,~] = cart2pol(xx,yy);
[th,r] = cart2pol(data(:,1),data(:,2));
% The hyperparameters are pretty reasonable, but here we fine tuned them,
% maximum likelihood style. 
% Not using the OOG points to get estimate smoothness, error.
%hyp2 = minimize(hyp, @gp, -100,@infExact, meanfunc, covfunc, likfunc, [data(~Loog,3),th(~Loog)],log10(r(~Loog)));
% Now fitting the data (inefficient, since a lot of the computation here is
% shared with the line above).
% Setting up a grid of points on which to evaluate the posterior
sqrtn = size(xx,1);
[x1,x2] = meshgrid(linspace(1,25,sqrtn),linspace(min(th),max(th),sqrtn));
x_star = [x1(:) x2(:)]; % column order: TF, theta
[m,s2,pm,ps2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, [data(:,3) th],log10(r),x_star);
%either error because x_star meshgrid can't be plotted, or error because zz
%th_star is too large. ???
%[m,s2,pm,ps2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, [data(:,3) th],log10(r),[zz(:) th_star(:)]);
ccs = [cos(x2(:,1)) sin(x2(:,1))];
[~,scalars] = gamutCheck([ccs zeros(size(ccs,1), 3-size(ccs,2))]', bkgndrgb, M, 'both');
gamutmask = repmat(log10(scalars'), 1, size(x2, 2));
predOOGpoints = m>gamutmask(:); % points that are predicted to be out of gamut
s2(predOOGpoints) = 0; % if the GP fit is outside of the gamut, don't even try to test there

s2 = reshape(s2,size(xx));
V = sqrt(xx.^2+yy.^2)-10.^sqrt(s2); % sqrt(s2) is the standard deviation of the log10 radial distance 
FV = isosurface(xx,yy,zz,V,-2.2);
h = patch(FV);
set(h,'FaceAlpha',.3,'EdgeColor','none')

%making the figure pretty
xlabel('L-Cone Contrast');
ylabel('M-Cone Contrast');
zlabel('Temporal Frequency (Hz)');
set(h,'FaceColor','yellow','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);
set(gcf,'Renderer','painters');
set(gca, 'zscale', 'log');
%set(gca, 'View', [135 12]);
set(gca, 'View', [225 12]);
axis vis3d;
drawnow;
end