%5/8/18 - Part of lmtfpaperfigures and I think also iterateandplotfiles. 
function lmtfModelClouds(subj, modelType, plotType)
try
    model_struct = load('C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Zack\IsoSamp\private\data\LMTF.mat');
catch
    model_struct = load('~/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat');
end
all_subj_models = getfield(model_struct, subj);
model_subset = getfield(all_subj_models, 'legacy');
specific_model = getfield(model_subset, modelType);
num_params = size(specific_model,1);
subplots = round(num_params/5,0);
subj_eccs = all_subj_models.eccs;
[phi,r] = cart2pol(abs(subj_eccs(:,1))/10,subj_eccs(:,2)/10);
figure;
if plotType == 1
    plot3(phi,r,specific_model(1,:), 'k.');
    axis vis3d;
    xlabel('Phi'); ylabel('R'); zlabel('Parameter Value');
    title('Xi_{Lum}');
    set(gca, 'View', [135 12]);
    figure;
    plot3(phi,r,specific_model(7,:), 'k.');
    axis vis3d;
    xlabel('Phi'); ylabel('R'); zlabel('Parameter Value');
    title('Xi_{RG}');
    set(gca, 'View', [135 12]);
    %     figure; hold on;
    %     for i = 1:num_params
    %         if i == 1 || i == 7
    %             continue;
    %         else
    %           %subplot(2,1,2); hold on;
    %           plot3(phi,r,specific_model(i,:), 'color', rand(1,3), 'linestyle', 'none','marker', '.');
    %         end
    %     end
    %     xlabel('Phi'); ylabel('R'); zlabel('Parameter Value');
    %     set(gca, 'View', [135 12]);
    %     axis vis3d;
elseif plotType == 2
    hold on;
    for i = 1:num_params
        subplot(5, subplots, i); hold on;
        [p_val,b] = plotmesh(specific_model(i,:));
        plot3(phi, r, specific_model(i,:), 'k.');
        if sum(isnan(b))>0
            title([subj, ',', int2str(i), '; regression coeff is nan']);
        elseif sum(isreal(b)) == 0
            title([subj, ',', int2str(i), '; regression coeff is imaginary']);
        else
            if p_val <= .05
                title_col = 'r';
            else
                title_col = 'k';
            end
            title([subj, ',', int2str(i), '; p-val: ', num2str(p_val,2)], 'color', title_col);
        end
    end
elseif plotType == 3
    hold on;
    plot3(phi, r, specific_model(1,:), 'k.');
    plotmesh(specific_model(1,:));
    xlabel('\Phi'); ylabel('R'); zlabel('Parameter Values'); title('Xi_{Lum}');
    set(gca, 'view', [225 12]);
    axis vis3d;
    set(gcf,'Renderer','painters');
    figure; hold on;
    plot3(phi, r, specific_model(7,:), 'k.');
    xlabel('\Phi'); ylabel('R'); zlabel('Parameter Values'); title('Xi_{RG}');
    set(gcf,'Renderer','painters');
    set(gca, 'view', [225 12]);
    axis vis3d;
    plotmesh(specific_model(7,:));
    
elseif plotType == 4
    % No longer just plotting the raw xi_lum and xi_coefficients
    % ploting the maximum predicted sensitivity.
    TF = 6.18; 
    xi_lum_coeffs = specific_model(1,:);
    xi_rg_coeffs = specific_model(7,:);
    [~,~,tmp] = tf_fiterr2([1; specific_model(2:6,1); 1; specific_model(8:13,1)],[1 1 TF 0]); % Height of fitted LUM TCSF at "TF" Hz assuming xi_lum = 1
    XiLUMScaleFactor = 1/tmp; % sensitivity peak, setting xi_lum to 1. Third output argument from ft_fiterr2 is predicted *threshold*.
    [~,~,tmp] = tf_fiterr2([1; specific_model(2:6,1); 1; specific_model(8:13,1)],[1 -1 TF 0]); % Height of fitted LUM TCSF at "TF" Hz assuming xi_lum = 1
    XiRGScaleFactor = 1/tmp; % sensitivity peak, setting xi_lum to 1. Third output argument from ft_fiterr2 is predicted *threshold*.

    % RG in r, phi
    hold on; set(gcf, 'Renderer', 'painters');
    plot3(phi, r, xi_rg_coeffs.*XiRGScaleFactor, 'k.','markersize', 20, 'markerfacecolor', 'k');
    [phimesh,rmesh] = meshgrid(-pi/2:pi/50:pi/2,0:.1:max(r));        
    [~, predprojsRG] = plotmesh(xi_rg_coeffs.*XiRGScaleFactor, phimesh, rmesh, 1);
    colormap(hot);
    color_min_3d = min(caxis);
    color_max_3d = max(caxis);
    xlabel('\phi');  ylabel('r');
    zlabel('Contrast sensitivity'); title('Chromatic contrast sensitivity');
    zlims = [min([xi_lum_coeffs.*XiLUMScaleFactor xi_rg_coeffs.*XiRGScaleFactor]*.9),... 
    max([xi_lum_coeffs.*XiLUMScaleFactor xi_rg_coeffs.*XiRGScaleFactor])*1.1];
    set(gca,'View',[133 15],'Zscale','log','zlim',zlims);
    
    % Luminance in r, phi
    figure; hold on;
    set(gcf,'Renderer','painters');
    plot3(phi, r, xi_lum_coeffs.*XiLUMScaleFactor, 'k.', 'markersize', 20, 'markerfacecolor', 'k');
    plotmesh(xi_lum_coeffs.*XiLUMScaleFactor, phimesh, rmesh, 1);
    xlabel('\phi'); ylabel('r');
    colormap(hot);
    caxis manual;
    caxis([color_min_3d color_max_3d]);
    zlabel('Contrast sensitivity');title('Luminance contrast sensitivity');
    set(gca,'View',[133 15],'Zscale','log','zlim',zlims);
    
    % RG in x,y
    [x,y] = meshgrid([0:.5:10], [-10:.5:10]);
    [phimesh, rmesh] = cart2pol(x,y);
    [~, predprojsRG] = plotmesh(xi_rg_coeffs.*XiRGScaleFactor, phimesh, rmesh, 0);
    color_min = min(predprojsRG);
    color_max = max(predprojsRG);
    figure; hold on; colormap(hot);
    h = surf(x,y,reshape(predprojsRG,size(phimesh)),ones(size(rmesh)));
    set(h,'EdgeColor','none', 'FaceAlpha', .6);
    axis image;
    caxis manual;
    caxis([color_min color_max]);
    set(gca,'View',[0 90], 'TickDir', 'out');
    set(gca, 'Xlim', [0 10], 'Ylim', [-10 10]);
    set(h,'CDataMode','auto')
    xlabel('Horizontal Position (DVA)');
    ylabel('Vertical Position (DVA)'); title('Xi_{RG}');
    contour(x,y,reshape(predprojsRG,size(phimesh)),log10([30 25 20 15 10 5]),'k-','linewidth',2)
    
    % Slopes
    %phi = linspace(-pi/2,pi/2,100);
    %[~, predprojsRG] = plotmesh(xi_rg_coeffs.*XiRGScaleFactor, [1;1]*phi, [1;2]*ones(1,length(phi)), 0);
    % Luminance in x, y
    figure; hold on; colormap(hot);
    [~, predprojsLUM] = plotmesh(xi_lum_coeffs.*XiLUMScaleFactor, phimesh, rmesh, 0);
    h = surf(x,y,reshape(predprojsLUM,size(phimesh)),ones(size(rmesh)));
    set(h,'EdgeColor','none', 'FaceAlpha', .6);
    axis image;
    caxis manual;
    caxis([color_min color_max]);
    set(gca, 'Xlim', [0 10], 'Ylim', [-10 10]);
    set(gca,'View',[0 90],'TickDir', 'out');
    set(h,'CDataMode','auto')
    xlabel('Horizontal Position (DVA)');
    ylabel('Vertical Position (DVA)'); title('Xi_{LUM}');
    contour(x,y,reshape(predprojsLUM,size(phimesh)),log10([30 25 20 15 10 5]),'k-','linewidth',2)
    % keyboard; set(gca,'Visible','off'); export_fig junk1 -png -transparent
    
    c = colorbar;
    c.TickDirection = 'out';
    tempTicks = [5 10 15 20 25 30];
    c.Ticks = log10(tempTicks);
    c.TickLabels = 10.^c.Ticks;
    drawnow;
    cdata = c.Face.Texture.CData;
    cdata(end,:) = uint8(.6 * cdata(end,:));
    c.Face.Texture.ColorType = 'truecoloralpha';
    c.Face.Texture.CData = cdata;
    drawnow;
    c.Face.ColorBinding = 'discrete';
else
    lsm = size(specific_model, 1);
    for i = 1:lsm
        coeff = specific_model(i, :);
        hold on; colormap(winter);
        set(gcf,'Renderer','painters');
        subplot(5,3,i);
        plot3(phi, r, coeff, 'k.', 'markersize', 20, 'markerfacecolor', 'k');
        [pval, ~] = plotmesh(coeff);
        xlabel('\phi'); ylabel('r');
        zlabel('l. c. s.'); title([int2str(i), ' ', num2str(pval)]);
        set(gca,'View',[133 15]);
        axis square;
    end
end

    function [p_val,predprojs] = plotmesh(param, phimesh, rmesh, makeplot)
        try
        [b,~,~,~,r_stats] = regress(log10(param'),[ones(numel(phi),1),r(:),r(:).*sin(2*phi(:)),r(:).*cos(2*phi(:))]);
        catch
            keyboard
        end
        predprojs = ([ones(size(rmesh(:))), rmesh(:),rmesh(:).*sin(2*phimesh(:)),rmesh(:).*cos(2*phimesh(:))]*b);
        if makeplot
            param_mesh = surface(phimesh,rmesh,reshape(10.^predprojs,size(rmesh)),'EdgeColor','none','FaceAlpha',.5);
        end
        p_val = r_stats(3);
    end
end