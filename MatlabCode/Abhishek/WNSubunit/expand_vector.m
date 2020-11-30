function reformed_vec = expand_vector(orig_vec,nrandnums_perchannel,mask,num_frames)
% This is a function that expands the 6 dimensional vector into the 300
% dimensional vector using mask
global nstixperside
subunits = nrandnums_perchannel;
basis = [];
n_elements = nstixperside^2;
reformed_vec = zeros(3*n_elements,num_frames);
for j = 1:num_frames
    for i = 1:subunits
        vec_r  = mask(1:n_elements);
        vec_g  = mask(n_elements+1:2*n_elements);
        vec_b  = mask(2*n_elements+1:3*n_elements);
        vec_r(vec_r ~= i) = 0; vec_r(vec_r == i) = orig_vec(i,j); % R
        vec_g(vec_g ~= subunits+1+i) = 0; vec_g(vec_g == subunits+1+i) = orig_vec(subunits+i,j); %G
        vec_b(vec_b ~= 2*(subunits+1)+i) = 0; vec_b(vec_b == 2*(subunits+1)+i) = orig_vec(2*subunits+i,j); % B
        vec = [vec_r;vec_g;vec_b];
        reformed_vec(:,j) = reformed_vec(:,j) + vec;
    end
end

end

