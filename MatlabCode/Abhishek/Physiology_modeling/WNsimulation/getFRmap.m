function non_lin = getFRmap(projST,projall,new_nbins)
n_spike = hist3(projST(:,1:2),{new_nbins,new_nbins});
n_raw = hist3(projall(:,1:2),{new_nbins,new_nbins});
n_spike = n_spike(2:end-1,2:end-1);
n_raw = n_raw(2:end-1,2:end-1);
n_raw(n_raw==0) = 1; % You can either enter 'NaN'or '1' to avoid division by zero
non_lin = n_spike./n_raw;
end

