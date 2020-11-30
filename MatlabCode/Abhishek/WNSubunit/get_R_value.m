function [R, model_fit] = get_R_value(non_lin, vec_est,vec,X,Y,hist_bins_cut)

% Extract the r value by comparing the model and the actual data
model_fit = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        model_fit(i,j) = vec_est(find(hist_bins_cut>(vec*[Y(i,j);X(i,j)]),1));
    end
end
R = corrcoef(non_lin(:),model_fit(:));
R = R(1,2);
end

