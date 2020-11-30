% to compute mutual information criterion as a function of mu and sigma

clear;
clc;

% min/max values for fmap
maxf = 250;
minf = -50;

% min/max values for sigma (std of fmap)
maxsigma = 100;
minsigma = 0.01;

frange = linspace(minf, maxf, 50);
sigrange = linspace(minsigma, maxsigma, 50);
[fxx, sigxx] = meshgrid(frange, sigrange);
sigxx(sigxx==0) = 0.001; % to avoid numerical problem, i.e. dividing by zero
fsigpairs = [fxx(:) sigxx(:)];
howmanypairs = length(fsigpairs);

%%
% min/max values for r (# of spikes)
% maxr = 50;
% minr = 0;
% rbins = minr:maxr;
% rbins(1) = 0.1;
% howmanyrbins = length(rbins);
% integralTable = zeros(howmanypairs, 1);
val = zeros(howmanypairs, 1);

for whichpair = 1: howmanypairs
    [whichpair howmanypairs]
    
    mu = fsigpairs(whichpair,1);
    sig = fsigpairs(whichpair,2);
    
    %
    fbins = linspace(mu-3.5*sig,mu+3.5*sig,100)';
    df = diff(fbins(1:2));
    maxf = max(fbins);
    maxr = max(20,maxf+3*sig+2*sqrt(abs(maxf+4*sig)));
    rbins = (0:maxr);
    [ff,rr] = meshgrid(fbins,rbins);
    
    [g,dg,ddg] = logexp1(ff);
    % ddLtrm1 = -bsxfun(@times,rbins,(ddg.*g-dg.^2)./g.^2);
    % ddLtrm2 = ddg;
    % ddL = ddLtrm1+ddLtrm2;
    %
    ddLtrm1 = -rr.*(ddg.*g-dg.^2)./g.^2;
    ddLtrm2 = ddg;
    ddL = ddLtrm1+ddLtrm2;
    pr = poisspdf(rr,g);
    pf  = normpdf(ff,mu,sig);
    
%     subplot(3,3,[2 3 5 6]);
%     imagesc(fbins,rbins,pr.*ddL.*pf);
%     axis xy;
%     title(sprintf('mu=%.2f, sig=%.2f', mu,sig));
%     
%     subplot(3,3,[1 4]);
%     plot(poisspdf(rbins,logexp1(mu)),rbins); axis tight;
%  ylabel('r');

%     subplot(3,3,8:9);
%     plot(fbins,normpdf(fbins,mu,sig)); axis tight;
%     xlabel('f');
%     
%     pause(0.05);
    
     val(whichpair) = sum(sum(pr.*ddL.*pf))*sig.^2*df;
    
%      val(whichpair) = sum(sum(pr.*ddL.*pf))*df;
    
    %%
    
    %     mubins1 = mu:sig/4:mu+3*sig;
    %     mubins2 = mu-sig:-sig/4:mu-3*sig;
    % %     mubins = linspace(
    %     mubins = [fliplr(mubins2)'; mubins1'];
    %     del_mu = mubins(2) - mubins(1);
    %     howmanymubins = length(mubins);
    %
    %     outerIntegral = zeros(howmanymubins, 1);
    %
    %     for whichmubin = 1:howmanymubins
    %         f = mubins(whichmubin);
    %         [g, dg, ddg] = logexp1(f);
    %
    %         innerIntegral = zeros(howmanyrbins, 1);
    %         for whichrbin = 1:howmanyrbins
    %             r = rbins(whichrbin);
    %             poissonLikeli = exp(-g)*g.^r./factorial(r);
    %             innerOtherTrm = -r*(ddg*g - dg.^2)./g.^2 + ddg;
    %             innerIntegral(whichrbin) = poissonLikeli*innerOtherTrm;
    %         end
    %
    % %         innerIntegral = dg.^2./g;
    %         normalprob = 1./(sqrt(2*pi)*sig).*exp(-0.5*(f-mu).^2./sig.^2);
    %         outerIntegral(whichmubin) = del_mu*normalprob*sum(innerIntegral);
    %     end
    %
    %     integralTable(whichpair) = sum(outerIntegral);
    
end

%%
% nsigs = length(sigrange);
% ss = fsigpairs(:,2);
% mm = fsigpairs(:,1);
% deltaruleSurf = ss.^2.*exp(mm + ss.^2./2);
% deltaruleSurf = reshape(deltaruleSurf,nsigs,[]);
% figure(2);
% surf(frange, sigrange,deltaruleSurf);
% xlabel('mu');
% ylabel('sigma'); 
% title('exp')

%%

% nsigs = length(sigrange);
% ss = fsigpairs(:,2);
% mm = fsigpairs(:,1);
% deltaruleSurf = ss.^2.*(exp(mm)./(exp(mm)+1)).^2./logexp1(mm);
% deltaruleSurf = reshape(deltaruleSurf,nsigs,[]);
% figure(2);
% surf(frange, sigrange, deltaruleSurf);
% xlabel('mu')
% ylabel('sigma'); title('us')


%%
% uncrtRule = reshape(ss.^2,nsigs,[]);
% figure(3);
% surf(frange, sigrange, uncrtRule);
% xlabel('mu')
% ylabel('sigma');


%%
% figure(2);surf(frange, sigrange, reshape(val, length(sigrange), []));
% xlabel('mu')
% ylabel('sigma');
% axis xy;

%%

% integral_in2D = reshape(integralTable, size(fxx,1), []);
musig2Dtable = [fsigpairs val];

datastruct2Dtable.frange = frange;
datastruct2Dtable.sigrange = sigrange;
datastruct2Dtable.fxx = fxx;
datastruct2Dtable.sigxx = sigxx;
datastruct2Dtable.val = reshape(val, length(sigrange), []);
save datastruct2Dtable datastruct2Dtable;

%%

load datastruct2Dtable.mat
figure(1); 
surf(datastruct2Dtable.frange, datastruct2Dtable.sigrange, reshape(datastruct2Dtable.val, length(datastruct2Dtable.sigrange), []));axis xy
xlabel('mu')
ylabel('sigma'); title('mi')


