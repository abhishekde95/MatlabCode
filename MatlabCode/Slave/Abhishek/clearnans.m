function out = clearnans(in)
col = size(in,2);
out = zeros(size(in));
for ii = 1:col
    data = in(:,ii);
    nan_idx = isnan(data);
    low_idx = find(nan_idx == 0,1);
    hi_idx = find(flipud(nan_idx) == 0,1);
    out(:,ii) = data;
end



end

