function plot_regression(xdata, ydata, x1, x2, yfit, xnew, ynew)
global gl
% x1 is in log10(TF) (from meshgrid)
% x2 is in radians (from meshgrid)
% yfit is in log10 (threhsold units)

xdata(:,1) = 10.^xdata(:,1);
ydata(:,1) = 10.^ydata(:,1);
xnew(:,1) = 10.^xnew(:,1);
ynew(:,1) = 10.^ynew(:,1);

figure(gl.HFIG); clf; axes;
figure_dims = [450 700];
screen_dims = get(0, 'ScreenSize');
set(gl.HFIG, 'position', [screen_dims(3)-figure_dims(1)-200 50 figure_dims]);

plot_helper(1, 10.^x1, x2, 10.^yfit, xdata, ydata, xnew, ynew);
plot_helper(2, 10.^x1, x2, 1./10.^yfit, xdata, 1./ydata, xnew, 1./ynew);

function plot_helper(k, x1, x2, yfit, xdata, ydata, xnew, ynew)
subplot(2,1,k); hold on;
% lay down the fitted surface
surf(x1, x2, reshape(yfit, size(x1)), 'edgecolor', 'none', 'linestyle', 'none');
xlabel('TF'); ylabel('\theta');
if k == 1
    zlabel('threshold');
else
    zlabel('sensitivity');
end
view(3); colormap(pmkmp([], 'linlhot')); shading interp;
set(gca, 'xscale', 'log', 'zscale', 'log');
set(gca, 'xlim', [min(xdata(:,1)) max(xdata(:,1))], 'ylim', [0 pi], ...
    'zlim', 10.^[floor(log10(min(ydata))) ceil(log10(max(ydata)))]);
set(gca, 'ytick', linspace(0, pi, 5), 'yticklabel', {'0' 'pi/4' 'pi/2', '3pi/4', 'pi'});
% dots represent measured data
plot3(xdata(:,1), xdata(:,2), ydata, 'o', 'markeredgecolor', 'k', 'markersize', 5, 'markerfacecolor', [.9 .9 .9]);
plot3(xnew(:,1), xnew(:,2), ynew, 'c*', 'markersize', 10);
