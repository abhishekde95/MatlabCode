%Part of IterateAndPlotFiles
function rawDataPlot(data, fpar, splitTitle, fit, zoom, symbol)
disp('raw data plot function');
%adjusting raw data to capture only TFs 0-10
if nargin < 6
    symbol = 'o';
    if nargin == 3
        fit = 0;
        zoom = 0;
    end
end
zoom_data = [];
if zoom
    for i = 1:length(data)
        if data(i,3) <= 5.200
            zoom_data = [zoom_data; data(i,:)];
        else
            continue
        end
    end
    data = zoom_data;
end

L = data(:,1);
M = data(:,2);
TF = data(:,3);
Loog = logical(data(:,4));

% Plotting the raw data
figure; axes; hold on;

plot3(L(~Loog),M(~Loog),TF(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5, 'Marker', symbol);
plot3(-L(~Loog),-M(~Loog),TF(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5, 'Marker', symbol);
plot3(L(Loog),M(Loog),TF(Loog),'ro','MarkerFaceColor','red','MarkerSize',5, 'Marker', symbol);
plot3(-L(Loog),-M(Loog),TF(Loog),'ro','MarkerFaceColor','red','MarkerSize',5, 'Marker', symbol);

% Adjusting the axes so we can see everything
if zoom
    lim = 0.14;
else
    lim = 0.97;
end

set(gca,'Xlim',1.1*[-lim lim]);
set(gca,'Zscale','log');
set(gca,'Ylim',1.1*[-lim lim]);
if zoom
    set(gca,'Zlim',[1 5]);
end
%axis square
xlabel('L-cone contrast');
ylabel('M-cone contrast');
zlabel('Temporal Frequency (Hz)');
set(gca,'View',[135 12]);
%axis square
axis vis3d

if fit
    % Plotting the fit
    if zoom
        [xx,yy,zz] = meshgrid(linspace(-.1,.1,30),...
            linspace(-.1,.1,30),...
            linspace(min(TF),5,30));  
    else
        [xx,yy,zz] = meshgrid(linspace(-1,1,50),...
            linspace(-1,1,50),...
            linspace(1,100,50));
    end
    
    % fits to one of the several rotations of the data that were tried.
    % f1 is always the filter in the direction of "besttheta". f2 is always
    % the filter in the orthogonal direction.
    % Order of model parameters
    xi_1 = fpar(1);
    zeta_1 = fpar(2);
    n1_1 = fpar(3);
    n2_1 = fpar(3)+fpar(4); % convention: n2 = n1+delta_n
    tau1_1 = 10^fpar(5);
    tau2_1 = 10^(fpar(5)+fpar(6)); % convention: tau2 = kappa*tau1
    xi_2 = fpar(7);
    zeta_2 = fpar(8);
    n1_2 = fpar(9);
    n2_2 = fpar(9)+fpar(10);
    tau1_2 = 10^fpar(11);
    tau2_2 = 10^(fpar(11)+fpar(12));
    theta = fpar(13);
    
    f1 = @(omega)xi_1*abs(((1i*2*pi*tau1_1.*omega+1).^-n1_1)-zeta_1*((1i*2*pi*tau2_1.*omega+1).^-n2_1));
    f2 = @(omega)xi_2*abs(((1i*2*pi*tau1_2.*omega+1).^-n1_2)-zeta_2*((1i*2*pi*tau2_2.*omega+1).^-n2_2));
    
    a = abs(f1(zz)); % luminance sensitivity
    b = abs(f2(zz)); % chromatic sensitivity
    mechs = [cos(theta) 1/sqrt(2); sin(theta) -1/sqrt(2)];
    V = zeros(size(xx));
    for i = 1:numel(xx)
        stimvector = [xx(i) yy(i)];
        V(i) = sqrt(sum((stimvector*mechs*diag([a(i) b(i)])).^2));
    end
    
    FV = isosurface(xx,yy,zz,V,1);
    h = patch(FV);
    %set(h,'FaceColor','green','EdgeColor','none');
    set(h,'FaceColor','green','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);
    set(gcf,'Renderer','painters');
    fv = tf_fiterr2(fpar,[L M TF],Loog); % error in fit
    title(num2str(fv));
end
% Good views
set(gca,'View',[135 12]);
%set(gca,'View',[225 12]);
lighting phong
set(gcf, 'name', splitTitle{1}, 'numbertitle', 'off');
drawnow;
end