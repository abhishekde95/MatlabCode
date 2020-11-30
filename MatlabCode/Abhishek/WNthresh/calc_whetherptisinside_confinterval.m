function [distance, rad, WwSTAsubdist] = calc_whetherptisinside_confinterval(STAsub,ptssub,vec,STCsub,nspikes,alpha)
% keyboard;
% [v,d] = eig(cov(ptssub-repmat(STAsub',size(ptssub,1),1)));
[v,d] = eig(STCsub);
Wwsub = sqrt(inv(d))*v';
Wwvec = Wwsub*(vec/norm(vec));
WwSTAsub = Wwsub*STAsub;
Wwvec_proj_STAsub = (dot(Wwvec,WwSTAsub)/(norm(Wwvec)^2)) * Wwvec;
% WwSTAsub'*(Wwvec/norm(Wwvec))
crit  = finv(1-alpha,3,nspikes-3);
crit = crit*((nspikes-1)*3)/(nspikes-3);
rad = sqrt(crit);
mdistance = sqrt((Wwvec_proj_STAsub-WwSTAsub)'*inv(nspikes*STCsub)*(Wwvec_proj_STAsub-WwSTAsub));
distance = norm(Wwvec_proj_STAsub - WwSTAsub);
Wwptssub = (Wwsub*ptssub')';
% tmprad = Wwptssub - repmat(WwSTAsub',size(Wwptssub,1),1);
% rad = max(sqrt(sum(tmprad.^2,2)));
WwSTAsubdist = norm(WwSTAsub);
% keyboard;
end

