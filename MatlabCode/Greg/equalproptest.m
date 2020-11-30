function [h,p] = equalproptest(counts,ns,alpha)
%   This function will perform a test
% to determine whether 'n' proportions
% are the same or different.
%
% function [h,p] = EqualPropTest(counts,ns,alpha)
%
% (See Snedecor and Cochran for details)

if (size(ns) ~= size(counts))
  error('The list of counts and ns have to have the same length');
end

if size(ns,2) == 1
    ns = ns';
    counts = counts';
end

props = counts./ns;

p = zeros(size(props,1),1);
h = zeros(size(props,1),1);

for i = 1:size(props,1)
    % From Snedecor and Cochran p.203
    pbar = sum(counts(i,:))/sum(ns(i,:));
    X2 = (sum(props(i,:).*counts(i,:)) - pbar*sum(counts(i,:)))/(pbar*(1-pbar));
    p(i) = 1-chi2cdf(X2,length(counts(i,:))-1);
    h(i) = p(i)<alpha;
end
% Testing
%ngroups = 5;
%n = 40;
% for i = 1:2000
%   counts = binornd(n,.2,ngroups,1);
%   [h,p] = EqualPropTest(counts,repmat(n,ngroups,1),.05);
%   ps = [ps p];
% end
% This is a meaningless change to the code to test out SVN
% And another