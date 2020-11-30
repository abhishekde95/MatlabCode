% For Analysing data stored in the structure gl (old file but works well)
 if ~exist('gl','var') 
     [fname, pathname] = uigetfile('*.mat', 'Select a data file','MultiSelect','on');
     if isequal(fname,0) || isequal(pathname,0)
         return
     end
    filename = [pathname fname];
    load(filename);
 end

global gl

new_idx = find_appro_index(gl);
prctcor = gl.cum_correct_count./(gl.cum_correct_count + gl.cum_incorrect_count);
prctcor = clearnans(prctcor);
new_x = logspace(log10(gl.stim.scale_lattice(1)),log10(gl.stim.scale_lattice(end)),101);
N_sf = numel(gl.stim.sf);
thresh = [];

for ii  = 1:N_sf
    sf_ind = [find(new_idx==N_sf*ii-1) find(new_idx==N_sf*ii)];
    
    % data and the psychometric function
    idx1 = ~isnan(prctcor(:,sf_ind(1)));
    [alpha1, beta1] = weibullFit(gl.stim.scale_lattice(idx1)', [gl.cum_correct_count(idx1,sf_ind(1)) gl.cum_incorrect_count(idx1,sf_ind(1))], 'mle');
    p1 = 1 +(0.5-1).*exp(-((gl.stim.scale_lattice(idx1)./alpha1).^beta1));
    q1 = 1 +(0.5-1).*exp(-((new_x./alpha1).^beta1));
    idx2 = ~isnan(prctcor(:,sf_ind(2)));
    [alpha2, beta2] = weibullFit(gl.stim.scale_lattice(idx2)', [gl.cum_correct_count(idx2,sf_ind(2)) gl.cum_incorrect_count(idx2,sf_ind(2))], 'mle');
    
    thresh = [thresh; alpha1 alpha2];
    q2 = 1 +(0.5-1).*exp(-((new_x./alpha2).^beta2));
    p2 = 1 +(0.5-1).*exp(-((gl.stim.scale_lattice(idx2)./alpha2).^beta2));
    
    % graphical representation of the threshold in the 2 color directions
    theta = pi/4:pi/2:(2*pi)-pi/4;
    rho = zeros(size(theta));
    if gl.color_ID{sf_ind(1)} == strcat('(L-M)-S_sf',num2str(ii))
        rho(1) = alpha2;
        rho(3) = alpha2;
        rho(2) = alpha1;
        rho(4) = alpha1;
    elseif gl.color_ID{sf_ind(2)} == strcat('(L-M)+S_sf',num2str(ii))
        rho(1) = alpha1;
        rho(3) = alpha1;
        rho(2) = alpha2;
        rho(4) = alpha2;
    end
    
    theta = [theta, theta(1)];
    rho = [rho, rho(1)];
    
    % plotting the figures
    fig_h = figure('Visible','off'); 
    
    subplot(121);plot(gl.stim.scale_lattice, prctcor(:,sf_ind(1)),'bo');hold on;
    plot(gl.stim.scale_lattice, prctcor(:,sf_ind(2)),'mo');
    plot(new_x', q1,'b');
    plot(new_x', q2,'m');
    xlabel('Contrast levels'), ylabel('Percentage correct choices');
    legend(gl.color_ID{sf_ind(1)}, gl.color_ID{sf_ind(2)},'location','SouthEast'); title('Psychometric Function');
    set(gca,'Xscale','log'); axis([gl.stim.scale_lattice(1) gl.stim.scale_lattice(end) 0 1]);
    
    subplot(122); 
    polar(theta,rho,'k-'); hold off;ylabel('S'),xlabel('L-M');
    title('Threshold plot');hold off;
    
    
    uicontrol(fig_h,'Position',[320 350 150 20],'Style','text',...
        'String','Spatial_frequency - ','Units','normalized');
    uicontrol(fig_h,'Position',[440 350 20 20],'Style','text',...
        'String',num2str(gl.stim.sf(ii)),'Units','normalized');
    uicontrol(fig_h,'Position',[330 50 80 20],'Style','text',...
        'String','Threshold 1 - ','Units','normalized');
    uicontrol(fig_h,'Position',[400 50 50 20],'Style','text',...
        'String',num2str(thresh(ii,1)),'Units','normalized');
    uicontrol(fig_h,'Position',[330 30 80 20],'Style','text',...
        'String','Threshold 2 - ','Units','normalized');
    uicontrol(fig_h,'Position',[400 30 50 20],'Style','text',...
        'String',num2str(thresh(ii,2)),'Units','normalized');
    fig_h.Visible = 'on';
end

