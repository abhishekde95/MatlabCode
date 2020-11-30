
% Set transfer function
gfun = @logexp1;
% gfun = @exp;

% Set to 1 if you want to generate a new table  
% (be sure to rename function below!)
generateTable = 1;

if generateTable
    % ===================
    nmus = 1000;
    nsigs = 1000;
    
    % Set range over mu values (mean of f)
    mumin = -100;
    mumax = 300;
    
    % Set range over sigma values  (stdev of f)
    sigmin = .01;
    sigmax = 200;
    
    
    % --- Compute 2D table of variance(mu,sigma) ---------
    
    % Make Gaussian pdf
    dx = .01;
    xx = (-5:dx:7)';
    pwts = normpdf(xx,0,1);
    pwts = pwts/sum(pwts);
    
    % Set grid over mu values
    dmu = (mumax-mumin)/(nmus-1);
    mus = (mumin:dmu:mumax)';
    
    % Set grid over sig (stdev) values)
    dsig = (log(sigmax)-log(sigmin))/(nsigs-1);
    sigs = exp(log(sigmin):dsig:log(sigmax))';
    
    % Compute first two moments across the grid
    M1 = zeros(nmus,nsigs);  % first moment
    M2 = zeros(nmus,nsigs);  % second moment
    for j = 1:nsigs
	mm = gfun(bsxfun(@plus,mus,xx'*sigs(j)));  % g evaluated for all mus, all x values, and one sigma
	M1(:,j) = mm*pwts;
	M2(:,j) = mm.^2*pwts;
    end
    VV_logexp1 = M2-M1.^2;
    
    %% Embed in a function, and save
    
    fvar_logexp1_lin = @(mu,sig)(interp2(mus',log(sigs),VV_logexp1',mu,log(sig),'*linear'));
    fvar_logexp1_nearest = @(mu,sig)(interp2(mus',log(sigs),VV_logexp1',mu,log(sig),'*nearest'));
    fvar_logexp1_cubic = @(mu,sig)(interp2(mus',log(sigs),VV_logexp1',mu,log(sig),'*cubic'));
    save fvar_logexp1_cubic fvar_logexp1_cubic;
    save fvar_logexp1_nearest fvar_logexp1_nearest;
    save fvar_logexp1_lin fvar_logexp1_lin;
    
    fmean_logexp1_lin = @(mu,sig)(interp2(mus',log(sigs),M1',mu,log(sig),'*linear'));
    fmean_logexp1_nearst = @(mu,sig)(interp2(mus',log(sigs),M1',mu,log(sig),'*nearest'));
    fmean_logexp1_cubic = @(mu,sig)(interp2(mus',log(sigs),M1',mu,log(sig),'*cubic'));
    save fmean_logexp1_cubic fmean_logexp1_cubic;
    save fmean_logexp1_nearst fmean_logexp1_nearst;
    save fmean_logexp1_lin fmean_logexp1_lin;
    
    % ==== END making table ========================
else
    % ==== load pre-made table =====================
    load fvar_logexp1_lin;
    load fvar_logexp1_cubic;
end

%% Compare accuracy: interpolated to quadrature method

mu =  20;
sig = 5;

% Evaluate integral numerically with 5000 grid points
dx = 10*sig/5000;
x = (mu-8*sig):dx:(mu+8*sig);
p = normpdf(x,mu,sig);
gg = gfun(x);
m1 = sum(gg.*p)*dx;
m2 = sum(gg.^2.*p)*dx;
vval = m2-m1.^2;

% Compare errors
% Err_lin_cub_quadrature = abs(vval-[fvar_logexp1_lin(mu,sig),...
%     fvar_logexp1_cubic(mu,sig),...
%     comptransformedPostVar_quadrature(gfun,mu,sig)])
Err_lin_cub_quadrature = abs(vval-[fvar_logexp1_lin(mu,sig), fvar_logexp1_nearest(mu,sig), fvar_logexp1_cubic(mu,sig)])

%% Compare Speed

% evaluate at some points
npt_tst = 400;
muvals = randn(npt_tst)*5+20;
sigvals = rand(npt_tst)*100;

tic;
fvar_logexp1_lin(muvals,sigvals);
toc;

tic;
fvar_logexp1_cubic(muvals,sigvals);
toc;

% tic;
% comptransformedPostVar_quadrature(gfun,muvals,sigvals);
% toc;

% ===================
% Take home message:
% ===================
% 1) Linear method is the fastest (by around 3x), but less accurate than cubic
% 2) Cubic method is 2nd fastest and *usually* most accurate
% 3) Quadrature method is slowest; is most accurate when function is smooth
%    (i.e., sigma << mu); least accurate when sigma>mu.  Doesn't have range
%    constraint the way the interpolated function does
