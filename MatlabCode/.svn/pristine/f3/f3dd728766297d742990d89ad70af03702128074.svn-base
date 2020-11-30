function [n_spike,n_raw,non_lin] = calcplanarspikestats(rawrgb_onbasisvec,stindices,new_nbins)
[n_spike,~] = hist3([rawrgb_onbasisvec(1,stindices)', rawrgb_onbasisvec(2,stindices)'],{new_nbins,new_nbins});
[n_raw,~] = hist3([rawrgb_onbasisvec(1,:)', rawrgb_onbasisvec(2,:)'],{new_nbins,new_nbins});
n_spike = n_spike(2:end-1,2:end-1);
n_raw = n_raw(2:end-1,2:end-1);
n_raw(n_raw==0) = 1; % You can either enter 'NaN'or '1' to avoid division by zero
non_lin = n_spike./n_raw;
blur = 1;
if blur
    non_lin = padarray(non_lin,[2 2],'replicate'); % pad the array with the border elements
    filt = fspecial('gaussian',5,1.0); % building a gaussian filter
    non_lin = conv2(non_lin,filt,'same'); % convolving a gaussian filter with the firing rate map
    non_lin = non_lin(3:end-2,3:end-2);
end
end

