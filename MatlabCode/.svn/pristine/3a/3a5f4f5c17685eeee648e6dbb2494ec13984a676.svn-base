function [plot_counter, U] = projontoSVD(spiketmpvec, rawtmpvec, plot_counter, figurename)

[U,S,V] = svd(spiketmpvec);
min_val = floor(min([U(:,1)'*rawtmpvec U(:,2)'*rawtmpvec])*100)/100;
max_val = ceil(max([U(:,1)'*rawtmpvec U(:,2)'*rawtmpvec])*100)/100;
num_bins = 11;
bin_interval = (max_val-min_val)/num_bins;
nbins1 = min_val:bin_interval:max_val;
nbins = linspace(prctile([U(:,1)'*rawtmpvec U(:,2)'*rawtmpvec],5), prctile([U(:,1)'*rawtmpvec U(:,2)'*rawtmpvec],95),numel(nbins1)+2);
xmin = min(nbins); xmax = max(nbins);
new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];

[n_spike,~] = hist3([spiketmpvec'*U(:,1) , spiketmpvec'*U(:,2)],{new_nbins,new_nbins});
[n_raw,~] = hist3([rawtmpvec'*U(:,1), rawtmpvec'*U(:,2)],{new_nbins,new_nbins});
blur = 1;
if blur
    filt = fspecial('gaussian',5,1.0); % building a gaussian filter
    n_spike = padarray(n_spike,[2 2],'replicate'); % pad the array with the border elements
    n_spike = conv2(n_spike,filt,'same'); % convolving a gaussian filter with the firing rate map
    n_raw = padarray(n_raw,[2 2],'replicate'); % pad the array with the border elements
    n_raw = conv2(n_raw,filt,'same'); % convolving a gaussian filter with the firing rate map
    n_raw(n_raw==0) = 1; % You can either enter 'NaN'or '1' to avoid division by zero
    non_lin = n_spike./n_raw;
    non_lin = non_lin(3:end-2,3:end-2);
end

[X,Y] = meshgrid(nbins,nbins);
figure(plot_counter), set(gcf,'Name',figurename); subplot(121); plot(U(:,1)'*rawtmpvec,U(:,2)'*rawtmpvec,'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
plot(U(:,1)'*spiketmpvec,U(:,2)'*spiketmpvec,'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); 
xlabel('SVD1'), ylabel('SVD2'); title('Proj onto the SVD vectors'); hold off;
subplot(122), contour(non_lin); set(gca,'XTick',[],'YTick',[]); xlabel('SVD2'), ylabel('SVD1'); title('Firing  Rate');
plot_counter = plot_counter + 1;
end

