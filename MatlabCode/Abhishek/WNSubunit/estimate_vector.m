function lik =  estimate_vector(input,projs,Lspike,nbins,str)
% Maximum likelihood estimate of the vector onto which the data stored in
% the 'projs vector' will be projected

eps =  1e-5;
a = input(3);
b = input(4);
c = input(5);

% input(1) - along S2 subunit cardinal direction
% input(2) - along S1 subunit cardinal direction
vec = [input(1) input(2)]; vec = vec/sqrt(vec*vec');
spike_proj = [projs(Lspike>0,1) projs(Lspike>0,2)]*vec'; % dot product of the spike triggered data onto the vector
raw_proj = projs*vec'; % dot product of the raw data onto the vector
ub = ceil(max([max(spike_proj) max(raw_proj)]));
lb = floor(min([min(spike_proj) min(raw_proj)]));
hist_bins = lb:mean(diff(nbins)):ub;
[X,Y] = meshgrid(nbins',nbins);
max_proj_val = vec(2)*X + vec(1)*Y;
ind = find(hist_bins > max(max_proj_val(:))+0.1,1);
hist_bins_cut = hist_bins(hist_bins>=-1 & hist_bins<=hist_bins(ind));
[spike_hist,~] = hist(spike_proj,hist_bins);
[raw_hist,~] = hist(raw_proj,hist_bins); raw_hist(raw_hist==0) = 1;

% maximum likelihood estimate
spike_hist_cut = spike_hist(hist_bins>=-1 & hist_bins<=hist_bins(ind));
raw_hist_cut = raw_hist(hist_bins>=-1 & hist_bins<=hist_bins(ind));

% Assuming that the relation between the firing rate and the projection value has an
% exponential relation
if (strcmp(str,'exp'))
    prob = a*exp(b*(hist_bins_cut+c));
elseif (strcmp(str,'quad'))
    prob = a*hist_bins_cut.^2 + b*hist_bins_cut + c;
end
prob(prob>=1-eps) = 1-eps;  % Ugly hack to avoid log10(0) = -Inf errors.
prob(prob<eps) = eps;

% Assuming the error values for a given x value is Bernoulli distributed
logP = log10([prob(:), 1-prob(:)]);
lik = spike_hist_cut'.*logP(:,1) + (raw_hist_cut-spike_hist_cut)'.*logP(:,2);
lik = (-1)*sum(lik);

end

