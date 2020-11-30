%5/8/18 part of lmtfpaperfigures. Never ended up being added to the paper.
function lmtfMushroomPlot(fulldata)
XY = unique(fulldata(:,[5 6]), 'rows');
Xs = (XY(:,1))';
Xs = mapminmax(Xs, 0, .9);
Ys = (XY(:,2)/100)';
Ys = mapminmax(Ys, .01, .9);
%plotting the actual data
figure('position', [100 100 1500 600]); hold on; axis off;
for i = 1:length(Ys)
    axes; set(gca,'position', [Xs(i), Ys(i) .1, .1]);
    eccs = fulldata(:,5) == XY(i,1) & fulldata(:,6) == XY(i,2);
    L = fulldata(eccs,1); M = fulldata(eccs,2); TF = fulldata(eccs,3);
    fpar = LMTFfitting(fulldata(eccs,:));
    [xx,yy,zz] = meshgrid(linspace(-max(abs(L)),max(abs(L)),35),...
        linspace(-max(abs(M)),max(abs(M)),35),...
        linspace(min(TF),max(TF),35));
    xi_1 = fpar(1); zeta_1 = fpar(2);
    n1_1 = fpar(3); n2_1 = fpar(3)+fpar(4); % convention: n2 = n1+delta_n
    tau1_1 = 10^fpar(5); tau2_1 = 10^(fpar(5)+fpar(6)); % convention: tau2 = kappa*tau1
    xi_2 = fpar(7); zeta_2 = fpar(8);
    n1_2 = fpar(9); n2_2 = fpar(9)+fpar(10);
    tau1_2 = 10^fpar(11);  tau2_2 = 10^(fpar(11)+fpar(12));
    theta = fpar(13);
    f1 = @(omega)xi_1*abs(((1i*2*pi*tau1_1.*omega+1).^-n1_1)-zeta_1*((1i*2*pi*tau2_1.*omega+1).^-n2_1));
    f2 = @(omega)xi_2*abs(((1i*2*pi*tau1_2.*omega+1).^-n1_2)-zeta_2*((1i*2*pi*tau2_2.*omega+1).^-n2_2));
    a = abs(f1(zz)).^-1;
    b = abs(f2(zz)).^-1;
    thtmp = atan2(yy,xx)-theta;
    rtmp = (b.*a)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    V = sqrt(xx.^2+yy.^2)-rtmp;
    FV = isosurface(xx,yy,zz,V,0);
    h = patch(FV);
    set(h,'FaceColor','green','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);
    set(gcf,'Renderer','painters');
    axis square; axis vis3d; 
    lim = max(max(abs([L M])));
    set(gca,'Xlim',1.1*[-lim lim]);
    set(gca,'Zscale','log');
    set(gca,'Ylim',1.1*[-lim lim]);
    set(gca,'View',[225 12]);
    axis off; drawnow;
end
end

% XY = [13 -35;13 35;20 -50;20 0;20 50;35 0;37 -92;37 92;50 -60;50 -50;50 -40;50 -20;50 0;50 20;50 40;50 50;50 60;80 0;100 -100;100 0;100 100];
% Xs = (XY(:,1))';
% Xs = mapminmax(Xs, 0, .8);
% Ys = (XY(:,2)/100)';
% [Ys, settings] = mapminmax(Ys, .01, .8);
% figure('position', [100 100 1500 600]); hold on; axis off;
% for i = 1:length(Ys)
%     axes;
%     set(gca,'position', [Xs(i), Ys(i) .2, .2]);
%     plot([0 1], [0 1], 'k');
%     axis vis3d;axis off;
% end