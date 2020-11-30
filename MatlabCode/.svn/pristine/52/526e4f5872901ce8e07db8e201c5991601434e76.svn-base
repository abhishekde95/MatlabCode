function [plot_counter,nbins,nbins1] = compute_projection_val_thresh(plot_counter,mode,subunits,nbins2,nbins3,flag,ch)

% Computes projection values on the basis vectors derived from the subunits
out = STPROJmod('return');
projs = out{1};
Lspike = out{2};
clear STPROJmod;
clear out;

if (flag == 1)
    labels = ['STA'; 'PC1'; 'PC2'; 'PC3'; 'PCi'];
    X_label = labels(ch(1),:);
    Y_label = labels(ch(2),:);
else
    X_label = 'S1 proj';
    Y_label = 'S2 proj';
end
if (mod(mode,3)==2)
    min_val = floor(min(projs(:))*100)/100;
    max_val = ceil(max(projs(:))*100)/100;
    num_bins = 15;
    bin_interval = (max_val-min_val)/num_bins; 
    nbins1 = min_val:bin_interval:max_val;
    % choose the bins such that they lie between 5 and 95 percentile
    nbins = linspace(prctile(projs(:),5), prctile(projs(:),95),numel(nbins1)+2);
    if (size(projs,2) == 2)
        figure(plot_counter); scatter(projs(:,1),projs(:,2)); hold on;
        scatter(projs(Lspike>0,1),projs(Lspike>0,2),'r'); 
        xlabel(X_label); ylabel(Y_label);       
        grid on;
        axis([min(projs(:,1))-0.1 max(projs(:,1))+0.1 min(projs(:,2))-0.1 max(projs(:,2))+0.1]);
        hold off;
        plot_counter = plot_counter + 1;
    end
else % in the current siuation if mode == 1
    nbins1 = nbins3;
    nbins = nbins2;
end
xmin = min(nbins); xmax = max(nbins);

if (size(projs,2) == 3)
    % using multivariate histogram if u have more than 2 subunits
    % Not perfect yet. Need to fix it.
    [count_raw,~,~,~] = histcn(projs,nbins,nbins,nbins);
    [count_st,~,mid_st,~] = histcn(projs(Lspike>0,:), nbins,nbins,nbins);
    count_raw(count_raw == 0) = 1;
    non_lin = count_st./count_raw;
    
    figure(plot_counter);
    subplot(1,subunits,1); Z = squeeze(sum(non_lin,3)); Z = Z/max(Z(:)); surfl(Z);
    xlabel('S2 proj'), ylabel('S1 proj'), zlabel('Non-linearity(S1-S2)');
    subplot(1,subunits,2); Z = squeeze(sum(non_lin,2)); Z = Z/max(Z(:)); surfl(Z);
    xlabel('S3 proj'), ylabel('S1 proj'), zlabel('Non-linearity(S1-S3)');
    subplot(1,subunits,3); Z = squeeze(sum(non_lin,1)); Z = Z/max(Z(:)); surfl(Z);
    xlabel('S3 proj'), ylabel('S2 proj'), zlabel('Non-linearity(S2-S3)');
    plot_counter = plot_counter + 1;
    
elseif (size(projs,2) == 2)
    % The program will run this code block if there are 2 subunits.
    % Converting the histogram into probability distribution and Estimating the non-linearity
    new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];
    [n_spike,out_spike] = hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{new_nbins,new_nbins});
    [n_raw,out_raw] = hist3([projs(:,1), projs(:,2)],{new_nbins,new_nbins});
    n_spike = n_spike(2:end-1,2:end-1);
    n_raw = n_raw(2:end-1,2:end-1);
    
    figure(plot_counter);
    % left hand side image for the spike trigerred ensemble image
    subplot(221); hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{nbins1,nbins1}); xlabel('S2 Proj'); ylabel('S1 Proj');
    title('Spike Triggered Ensemble'); set(gcf,'renderer','opengl'); set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    subplot(222); hist3([projs(:,1), projs(:,2)],{nbins1,nbins1}); xlabel(X_label); ylabel(Y_label);
    title('Raw Ensemble'); set(gcf,'renderer','opengl'); set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    subplot(224); imagesc([min(new_nbins) max(new_nbins)],[min(new_nbins) max(new_nbins)], n_raw); xlabel(X_label); ylabel(Y_label); title('Raw Ensemble');
    subplot(223); imagesc([min(new_nbins) max(new_nbins)],[min(new_nbins) max(new_nbins)], n_spike); xlabel(X_label); ylabel(Y_label); title('Spike Triggered Ensemble');
    plot_counter = plot_counter + 1;
    n_raw(n_raw==0) = 1; % You can either enter 'NaN'or '1' to avoid division by zero
    non_lin = n_spike./n_raw;
    non_lin_pa = padarray(non_lin,[2 2],'replicate'); % pad the array with the border elements
    filt = fspecial('gaussian',5,1.0); % building a gaussian filter
    non_lin_blurred = conv2(non_lin_pa,filt,'same'); % convolving a gaussian filter with the firing rate map
    non_lin_blurred = non_lin_blurred(3:end-2,3:end-2);
        
    figure(plot_counter);
    subplot(2,2,1), surfl(nbins',nbins,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
    xlabel(Y_label), ylabel(X_label); title('Firing  Rate');
    subplot(2,2,2), plot(nbins,sum(non_lin,1)/max(sum(non_lin,1)),'LineWidth',2);   % plotting the normalised marginal non-linearity
    xlabel(Y_label), ylabel('Firing Rate'); title('Projection onto S2 vector');
    axis([nbins(1) nbins(end) 0 1]);
    subplot(2,2,3),plot(nbins,sum(non_lin,2)/max(sum(non_lin,2)),'LineWidth',2);
    ylabel('Firing Rate'), xlabel(X_label); title('Projection onto S1 vector');
    axis([nbins(1) nbins(end) 0 1]); %view([-90 90]);
    subplot(2,2,4), imagesc([xmin xmax],[xmin xmax],non_lin); 
    xlabel(Y_label), ylabel(X_label); title('2-D Firing  Rate Intensity Plot');
    plot_counter = plot_counter + 1;
    
    figure(plot_counter);
    subplot(231), imagesc([xmin xmax],[xmin xmax],non_lin);
    xlabel(Y_label), ylabel(X_label); title('FR w/o the filter');
    subplot(232), imagesc([xmin xmax],[xmin xmax],non_lin_blurred); 
    xlabel(Y_label), ylabel(X_label); title('FR after the filter');
    subplot(234), surfl(nbins',nbins,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
    xlabel(Y_label), ylabel(X_label); 
    subplot(235), surfl(nbins',nbins,non_lin_blurred); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin_blurred(:)) max(non_lin_blurred(:))]);
    xlabel(Y_label), ylabel(X_label);
    subplot(233),contour(nbins',nbins,non_lin_blurred);set(gca,'YDir','Reverse');
    xlabel(Y_label), ylabel(X_label); title('Contour plot');
    plot_counter = plot_counter + 1;
    
elseif (size(projs,2) == 1)
    
    % The program will run this code block if there is only 1 subunit.
    figure(plot_counter);
    % left hand side image for the spike trigerred ensemble image
    subplot(1,2,1); hist(projs(Lspike>0),nbins1);
    hold on;
    xlabel('Proj values'); ylabel('frequency'); title('Spike Triggered Ensemble'); hold off;
    subplot(1,2,2);hist(projs ,nbins1)
    hold on;
    xlabel('Proj values'); ylabel('frequency'); title('Raw Ensemble'); hold off;
    plot_counter = plot_counter + 1;
    
    % Converting the histogram into probability distribution and Estimating the non-linearity
    [n_spike,out_spike] = hist(projs(Lspike>0),nbins);
    [n_raw,out_raw] = hist(projs,nbins);
    n_raw(n_raw==0) = 1; % to avoid division by zero
    non_lin = n_spike./n_raw;
    
    figure(plot_counter);
    plot(nbins,sum(non_lin,1)/max(sum(non_lin,1)),'LineWidth',2);   % plotting the normalised marginal non-linearity
    xlabel('Proj values'), ylabel('Firing Rate');
    axis([nbins(1) nbins(end) 0 1]);
    plot_counter = plot_counter + 1;
end


end

