function initargs = define_basis_vec(nrandnums_perchannel,stro,basis_vec,new_vec,flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global nframesidx maxT;
if (flag == 1)
    nframestotal = sum(stro.trial(:,nframesidx));
    vec1 = reshape(new_vec(:,1),[nrandnums_perchannel*3 maxT]);
    vec1_rev = flipdim(vec1,2);
    vec2 = reshape(new_vec(:,2),[nrandnums_perchannel*3 maxT]);
    vec2_rev = flipdim(vec2,2);
    initargs = {[vec1_rev(:), vec2_rev(:)] , 0, nframestotal, [nrandnums_perchannel 3 maxT]};
    clear vec1 vec2 vec1_rev vec2_rev;
else
    
    subunits = nrandnums_perchannel;
    basis = [];
    frame_len = size(basis_vec,2);
    for i = 1:subunits
        vec = zeros(size(basis_vec));
        vec(i,:) = basis_vec(i,:);
        vec(subunits + i,:) = basis_vec(subunits + i,:);
        vec(2*subunits + i,:) = basis_vec(2*subunits + i,:);
        % reversing the vector, flipping the vector across columns
        vec_rev = flipdim(vec,2);
        tmp =  vec_rev(:);
        tmp = tmp/sqrt(sum(tmp.*tmp)); % normalised basis vectors
        basis = [basis, tmp];
        clear tmp vec_rev vec;
    end
    nframestotal = sum(stro.trial(:,nframesidx));
    initargs = {basis , 0, nframestotal, [nrandnums_perchannel 3 maxT]};
end


end
