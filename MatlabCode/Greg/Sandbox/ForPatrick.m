% PLOTTING A REGULAR TETRAHEDRON

% Equilateral triangle

a = [1 0];
b = [cosd(120) sind(120)];
c = [cosd(240) sind(240)];

% Dotproducts are maximally negative and norms are all '1'
[a; b; c] * [a;b;c]'

% Triangle is centered at the origin
mean([a; b; c])

% % Plotting
% figure; axes; hold on;
% plot([a(1) b(1)],[a(2) b(2)],'k-');exit
% plot([a(1) c(1)],[a(2) c(2)],'k-');
% plot([b(1) c(1)],[b(2) c(2)],'k-');
% axis equal;
% 

% Now the tetrhedron
% Here's what I want the outproduct matrix to look like:
A = eye(4);
A(A==0) = -1/3;

B = chol(A);
a = B([1 2 3],1);
b = B([1 2 3],2);
c = B([1 2 3],3);
d = B([1 2 3],4);

% It's centered at the origin
mean([a';b';c';d']);

% Plotting the tetrahedron
figure; axes; hold on;
plot3([a(1) b(1)],[a(2) b(2)],[a(3) b(3)],'k-');
plot3([a(1) c(1)],[a(2) c(2)],[a(3) c(3)],'k-');
plot3([a(1) d(1)],[a(2) d(2)],[a(3) d(3)],'k-');
plot3([b(1) c(1)],[b(2) c(2)],[b(3) c(3)],'k-');
plot3([b(1) d(1)],[b(2) d(2)],[b(3) d(3)],'k-');
plot3([c(1) d(1)],[c(2) d(2)],[c(3) d(3)],'k-');

set(gca,'XLim',[-1 1],'YLim',[-1 1],'Zlim',[-1 1]);
axis vis3d;

% Finding the midpoints of the legs.
mp = [];
for i = 1:6
    switch(i)
        case (1)
            mp(i,:) = mean([a';b']);
        case (2)
            mp(i,:) = mean([a';c']);
        case (3)
            mp(i,:) = mean([a';d']);
        case (4)
            mp(i,:) = mean([b';c']);
        case (5)
            mp(i,:) = mean([b';d']);
        case (6)
            mp(i,:) = mean([c';d']);

    end
    plot3(mp(i,1),mp(i,2),mp(i,3),'m*');
end

% Connecting the midpoints of the legs
for i = 1:4
    L = mp*B([1 2 3],i) > 0;
    tmp = mp(L,:);
    tmp = [tmp; tmp(1,:)];  % wrap
    plot3(tmp(:,1),tmp(:,2),tmp(:,3),'g-');
end

% Finding the center of each face
for i = 1:4
   tmp = B([1 2 3],:);
   tmp(:,i) = [];
   fp = mean(tmp,2);
   plot3(fp(1),fp(2),fp(3),'k*');
end


% The origin
 plot3(0,0,0,'ro')
 
%%
% PLAYING WITH DELAUNAY TRIANGULATION/VORONOI DIAGRAMS
 
x = unifrnd(0,1,10,1);
y = unifrnd(0,1,10,1);
dt = DelaunayTri(x,y)
triplot(dt);
[V,C] = voronoiDiagram(dt)


x = unifrnd(0,1,10,1);
y = unifrnd(0,1,10,1);
[vx,vy] = voronoi(x,y)
h = plot(vx,vy,'-',x,y,'.')

% -----------
% Simulation
% -----------
x = unifrnd(0,1,5,1);
y = unifrnd(0,1,5,1); 
x = [1 1 0 0]';
y = [1 0 0 1]';
th = unifrnd(0,2*pi);
trans = [cos(th) -sin(th); sin(th) cos(th)];
trans = normrnd(0,3,2,2);

%neuron = inline('abs(2*(x-.5).^2) + y','x','y');
neuron = inline('1./(1+exp(-1*((x-.5)+(y-.5))))','x','y');
%neuron = inline('(1./(1+exp(-5*((x-.5)+2*(y-.5)))))','x','y');
%neuron = inline('x+y','x','y');

[tmpx, tmpy] = meshgrid([0:.01:1],[0:.01:1]);
tmp_in = [tmpx(:) tmpy(:)];
tmp_out = trans*tmp_in';
tmpx = reshape(tmp_out(1,:),size(tmpx));
tmpy = reshape(tmp_out(2,:),size(tmpy));

figure; axes; hold on;
fr = neuron(tmpx,tmpy);
colormap(jet(64));
image([0:.01:1], [0:.01:1], fr./max(fr(:))*64);
set(gca,'YLim',[0 1],'XLim',[0 1]);

for i = 1:100
    xy = [x,y]*trans';
    
    plot(x,y,'k*')
    fr = poissrnd(neuron(xy(:,1),xy(:,2)));
    
    dt = DelaunayTri(x,y);
    [V,C] = voronoiDiagram(dt);
    v = unique(V(2:end,:),'rows');
    Linvalid = v(:,1) < 0 | v(:,2) < 0 | v(:,1) > 1 | v(:,2) > 1;
    v(Linvalid,:) = [];
    tmp = [];
    for j = 1:size(v,1)
        dists = [x y]-repmat(v(j,:),size(x,1),1);
        dists = sqrt(sum(dists.^2,2));
        L = softEq(dists, min(dists));
        if (sum(L) == 1)
            keyboard
        end
        tmp = [tmp; var(fr(L))./mean(fr(L))];
    end
    idx = find(tmp == max(tmp),1);
    
    x = [x; v(idx,1)];
    y = [y; v(idx,2)];
    drawnow;
    pause(.1)
end

% 
% % 3D
% x = unifrnd(0,1,10,3);
% dt = DelaunayTri(x(:,1),x(:,2),x(:,3));
% [V,C] = voronoiDiagram(dt);
% figure; axes; hold on;
% plot3(x(:,1),x(:,2),x(:,3),'k*');
% plot3(V(2:end,1),V(2:end,2),V(2:end,3),'m*');
% set(gca,'Xlim',[0 1],'Ylim',[0 1],'Zlim',[0 1])\
%%
% Finding the intersection of 2 lines
x1 = 4;
y1 = 6;
x2 = 8;
y2 = -10;

x3 = 9;
y3 = 1;
x4 = -1;
y4 = 2;

m1 = (y2-y1)/(x2-x1);
b1 = y1-m1*x1
m2 = (y4-y3)/(x4-x3);
b2 = y3-m2*x3


tmp = linspace(-1,10,30);
figure; axes; hold on;
plot(tmp, tmp*m1+b1)
plot(x1,y1,'m*')
plot(x2,y2,'m*')
plot(tmp, tmp*m2+b2)
plot(x3,y3,'b*')
plot(x4,y4,'b*')

xint = (b2-b1)/(m1-m2);
yint = m1*xint+b1;
plot(xint,yint,'g*');


%%
% Trying again.
% Bootstrap to identify regions of the stimulus space that are
% poorly constrained (by some criterion).

neuron = inline('10./(1+exp(-1*(-10*(x+.5)-3*y)))','x','y');
%neuron = inline('max(0,100*(x-.5))','x','y');
%neuron = inline('(10./(1+exp(-5*((x-.5)+2*(y-.5).^2))))','x','y')
noisesigma = 1;

% Plotting the firing rate function
figure(1); clf; axes; hold on;
[tmpx, tmpy] = meshgrid([-1:.01:1],[-1:.01:1]);
fr = neuron(tmpx,tmpy);
colormap(jet(64));
image([-1:.01:1], [-1:.01:1],(fr-min(fr(:)))./(max(fr(:))-min(fr(:)))*64);
set(gca,'YLim',[-1 1],'XLim',[-1 1]);


xy = fullfact ([3 3])-2;
x = xy(:,1)';
y = xy(:,2)';
x = [-1 1 1 -1 0];
y = [1 1 -1 -1 0];

%x = .25*cos(linspace(0,2*pi-pi/4,8));
%y = .25*sin(linspace(0,2*pi-pi/4,8));
fr = neuron(x,y)+normrnd(0,noisesigma,1,length(x));

smoothingparameter = .5;
niter = 40;
for i = 1:niter
    
    st = tpaps([x;y], fr,smoothingparameter);
    
    figure(1);
    plot(x,y,'w.');
    
    % Bootstrapping residuals
    %     preds = fnval(st,[x;y]);
    %     resid = preds-fr;
    %     nboot = 200;
    %     for i = 1:nboot
    %         bootfr = preds+resid(unidrnd(length(resid),length(resid),1));
    %         bootst(i) = tpaps([x;y], bootfr);
    %     end
    
    % Bootstrapping raw data
    
    nboot = 100;
    clear bootst;
    j = 1;
    while j < nboot
        bootidxs = unidrnd(length(x),length(x),1);
        try
            bootst(j) = tpaps([x(bootidxs);y(bootidxs)], fr(bootidxs),smoothingparameter);
            j = j +1;
        catch
            % Do nothing
        end
    end
    
    
    %  Plotting the stack of fits
    %  figure; axes; hold on;
    %  for i = 1:length(bootst)
    %      fnplt(bootst(i));
    %  end
    
    % Getting a map of errors in the predictions
    [tmpx, tmpy] = meshgrid([-1:.05:1],[-1:.05:1]);
    for j = 1:length(bootst)
        gridpreds(j,:) = fnval(bootst(j),[tmpx(:)';tmpy(:)']);
        err(j,:) = gridpreds(j,:) - fnval(st,[tmpx(:)';tmpy(:)']);
    end
    
    figure(2);
    subplot(2,2,1);
    mesh(reshape(sum(err).^2,size(tmpx)))
    subplot(2,2,2);
    mesh(reshape(sum(err.^2),size(tmpx)))
    subplot(2,2,3);
    tstat = abs(mean(err)./std(err));
    mesh(reshape(tstat,size(tmpx)))
    drawnow;
    
    % Picking a new place to go
    %   idx = find(sum(err).^2 == max(sum(err).^2));
    idx = find(tstat == max(tstat));
    x = [x, tmpx(idx)];
    y = [y, tmpy(idx)];
    fr = [fr, neuron(tmpx(idx),tmpy(idx))+normrnd(0,noisesigma)];
end

figure; axes; hold on;
plot3(x,y,fr,'ko','MarkerFaceColor','black')
% Plotting the true firing rate function
[tmpx, tmpy] = meshgrid([-1:.01:1],[-1:.01:1]);
mesh(tmpx,tmpy,neuron(tmpx,tmpy),ones(size(tmpy)));
st = tpaps([x;y], fr);
fnplt(st);

%%
% Trying again.
% Leave one out and fit models.   Create a "surprise" function which is the
% difference between each data point and the leave-one-out prediction at
% that location.  Go to the peak of the surprise function and find a local
% minimum in the proximity surface near it.
options = optimset('MaxIter',1000,'Tolx',10^-10,'TolFun',10^-10,'MaxFunEvals',1000);
%neuron = inline('5*(x+1)+2*(y+1)','x','y');

neuron = inline('max(0,10*(x-.5)+30*(y+.5))','x','y');
%neuron = inline('60*x.^2+60*y.^2+10.*x.*y','x','y');

%neuron = inline('(50./(1+exp(-5*((x-.5)+2*(y-.5).^2))))','x','y');
%neuron = inline('100./(1+exp(-1*(+5*(x-.5)-3*y)))','x','y');
%neuron = inline('10*max(0,-(x-.5)+2*y).^2','x','y');

niter = 100;  % number of stimuli
noisesigma = -1;
USESURPRISE = 0;
INITIALSMOOTHINGPARAMETER = .1;
DELTASMOOTHING = .01;
nreps = 5;  % Number of simulated experiments
PLOTFIGS = 1;
UPDATEDISPLAY = 1;
DWTHRESH = 2;
MSEs = [];
GRIDINC = .1;
if (PLOTFIGS)
    RESIDFIG = figure;
end
if (UPDATEDISPLAY)
    SAMPLINGFIG = figure;
end
for rep = 1:nreps
    disp('--------')
    disp(['Simulation ',num2str(rep)]);
    
    smoothingparameter = INITIALSMOOTHINGPARAMETER;  % Need to reset it
    % Plotting the firing rate function
    
    if (UPDATEDISPLAY)
        figure(SAMPLINGFIG); clf; axes; hold on;
        [tmpx, tmpy] = meshgrid([-1:.01:1],[-1:.01:1]);
        fr = neuron(tmpx,tmpy);
        colormap(jet(64));
        image([-1:.01:1], [-1:.01:1],(fr-min(fr(:)))./(max(fr(:))-min(fr(:)))*64);
        set(gca,'YLim',[-1 1],'XLim',[-1 1]);
    end
    
    %xy = fullfact ([3 3])-2;
    %x = xy(:,1)';
    %y = xy(:,2)';
    x = [-1 1 1 -1];
    y = [1 1 -1 -1];
    
    if (noisesigma>=0)
        fr = neuron(x,y)+normrnd(0,noisesigma,1,length(x));
    else
        fr = poissrnd(max(0,neuron(x,y)));
    end
    
    
    for i = 1:niter
        % Computing the "surprise" function
        st = tpaps([x;y], fr,smoothingparameter); 
        surprise = zeros(length(x),1);
        
        if (USESURPRISE == 1 | USESURPRISE == 0)
            for j = 1:length(x)
                tmpx = x; tmpy = y; tmpfr = fr;
                tmpx(j) = []; tmpy(j) = []; tmpfr(j) = [];
                fit = tpaps([tmpx;tmpy], tmpfr,smoothingparameter);
                fullmodelpred = fnval(st,[x(j);y(j)]);
                reducedmodelpred = fnval(fit,[x(j);y(j)]);
                surprise(j) = abs(max(0,fullmodelpred)-max(0,reducedmodelpred));
              
   
                %    surprise(j) = sqrt(abs(max(0,fullmodelpred)-max(0,reducedmodelpred)));
                
                %    surprise(j) = surprise(j)-(b(1)*max(0,fullmodelpred)+b(2));  %
                %    didn't work
            end
            if (USESURPRISE == 0)
                surprise = surprise(randperm(length(surprise)));
            end
            % Splining the surprise function
            surprisespline = tpaps([x;y], surprise');
        else % USESURPRISE == 2
            nboot = 100;
            clear bootst;
            j = 1;
            while j < nboot
                bootidxs = unidrnd(length(x),length(x),1);
                try
                    bootst(j) = tpaps([x(bootidxs);y(bootidxs)], fr(bootidxs),smoothingparameter);
                    j = j +1;
                catch
                    % Do nothing
                end
            end
            [tmpx, tmpy] = meshgrid([-1:GRIDINC:1],[-1:GRIDINC:1]);
            for j = 1:length(bootst)
                err(j,:) = fnval(bootst(j),[tmpx(:)';tmpy(:)']) - fnval(st,[tmpx(:)';tmpy(:)']);
            end
          %  surprisegrid = reshape(abs(mean(err)./std(err)),size(tmpx));
          
            surprisegrid = reshape(sum(sqrt(abs(err))),size(tmpx));
            surprisespline = tpaps([tmpx(:), tmpy(:)]', surprisegrid(:)');
        end
        
        
        %  if (PLOTFIGS)
        %      figure; axes; hold on; fnplt(surprisespline);
      %      plot3(x,y,surprise','k.');
      %  end
      
      % Adjusting smoothing parameter if need be
      [~,rtidx] = sort(fnval(st,[x;y]));
      resid = fr-max(0,fnval(st,[x;y]));
      %[h,p] = runstest(resid)
      dw = sum((resid(2:end)-resid(1:end-1)).^2)./sum(resid.^2);
      if (dw > DWTHRESH & length(x) > 10)
          smoothingparameter = smoothingparameter+(1-smoothingparameter)*DELTASMOOTHING;
      end
      
      
      if (PLOTFIGS)
          figure(RESIDFIG);
          subplot(2,2,1);
          plot(fnval(st,[x;y]),resid,'k.');
       
          title(['dw = ',num2str(dw,2),' Smoothing: ',num2str(smoothingparameter,3)]);
          ylabel('residual'); xlabel('predicted value');
          subplot(2,2,2); cla; hold on;
          plot(fnval(st,[x;y]), abs(resid),'k.');
          b = [fnval(st,[x;y])' ones(length(resid),1)]\abs(resid)';
          plot([0 max(fnval(st,[x;y]))],[b(2) b(1)*max(fnval(st,[x;y]))+b(2)]);
          ylabel('squared residual'); xlabel('predicted value');
          
          subplot(2,2,3);
          fnplt(surprisespline);
          set(gca,'View',[0 90]);
          title('surprise');
          subplot(2,2,4);
          residspline = tpaps([x;y],fr-max(0,fnval(st,[x;y])));
          fnplt(residspline);
          set(gca,'View',[0 90]);
          title('residuals');
      end
      
      tmpgrid = [-1:GRIDINC:1];
      
      surprisegrid = [];
      for j = 1:length(tmpgrid)
          for k = 1:length(tmpgrid)
              surprisegrid(j,k) = fnval(surprisespline,[tmpgrid(j); tmpgrid(k)]);
          end
      end

%         
        [row,col] = ind2sub(size(surprisegrid),find(surprisegrid(:) == max(surprisegrid(:))));
        mspx = tmpgrid(row)*.99+normrnd(0,.0001);
        mspy = tmpgrid(col)*.99+normrnd(0,.0001);
        

%         if (USESURPRISE)
%             msp = find(abs(surprise) == max(abs(surprise)));
%             mspx = x(msp)*.99+normrnd(0,.0001);
%             mspy = y(msp)*.99+normrnd(0,.0001);
%         else
%             mspx = unifrnd(-1,1);
%             mspy = unifrnd(-1,1);
%         end
        % This function minimization is still having some problems.
        
        %disp(['starting at :',num2str([mspx mspy])]);
      %  [point,fval,exitflag,output] = fminsearch(@(startpoint) proximity(startpoint,x,y),[mspx mspy],options); 
        [point,fval,exitflag,output] = fminsearch('proximity',[mspx mspy],options,x, y); 
        %point = [mspx mspy];
        
        if (UPDATEDISPLAY)
            figure(SAMPLINGFIG);
            plot(x,y,'ko','MarkerSize',5,'MarkerFaceColor','black');
        end       
        
        
        drawnow;
        %pause;
        
        x = [x, point(1)];
        y = [y, point(2)];
        
        if (noisesigma>=0)
            fr = [fr, neuron(point(1),point(2))+normrnd(0,noisesigma)];
        else
            fr = [fr, poissrnd(neuron(point(1),point(2)))];
        end
    end  % End of main loop
    
    % Evaluating the fit
    [tmpx, tmpy] = meshgrid([-1:.05:1],[-1:.05:1]);

%    st = tpaps([x;y], fr);
%    predvals = zeros(size(tmpx));
%    for i = 1:length(tmpx(:))
%        predvals(i) = max(0,fnval(st,[tmpx(i);tmpy(i)]));
%    end
    % Evaluating the quality of the fit on the basis of a polynomial
    % regression fit (assuming Poisson errors?)
%    designmat = [x; y; x.*y; x.^2; y.^2 ]';
 %   designmat = [x; y; ]';
    
%    bglm = regress(fr',[ones(length(x),1) designmat]);
%    bglm = glmfit(designmat ,fr','poisson','link','identity')
%    predvals = [ones(length(tmpx(:)),1) tmpx(:) tmpy(:) tmpx(:).*tmpy(:) tmpx(:).^2 tmpy(:).^2]*bglm;
%    predvals = [ones(length(tmpx(:)),1) tmpx(:) tmpy(:)]*bglm;  
    
    % Spline fit to data for assesing MSE
    st = tpaps([x;y], fr);
    predvals = zeros(size(tmpx));
    for i = 1:length(tmpx(:))
        predvals(i) = max(0,fnval(st,[tmpx(i);tmpy(i)]));
    end
    
    truth = reshape(neuron(tmpx(:),tmpy(:)),size(tmpx));
    MSE = mean((truth(:)-predvals(:)).^2);
    MSEs(rep) = MSE;
    if (PLOTFIGS)
        figure; axes; hold on;
        plot3(x,y,fr,'ko','MarkerFaceColor','black')
        % Plotting the true firing rate function
        mesh(tmpx,tmpy,truth,ones(size(tmpy)));
        % predvals = [ones(length(tmpx(:)),1) tmpx(:) tmpy(:) tmpx(:).*tmpy(:) tmpx(:).^2 tmpy(:).^2]*bglm;
        % predvals = [ones(length(tmpx(:)),1) tmpx(:) tmpy(:)]*bglm;
        mesh(tmpx,tmpy,predvals,3*ones(size(tmpy)));

        
        title(['MSE = ',num2str(MSE)]);
    end
    
    % same number of points, but uniform sampling
end % reps loop

if (nreps > 1)
    figure; axes; hold on;
    hist(MSEs);
    plot(mean(MSEs)+[-std(MSEs)/sqrt(nreps) std(MSEs)/sqrt(nreps)],[0 0],'g-','Linewidth',2);
    plot(mean(MSEs),0,'m*')
    title(['Use surprise: ',num2str(USESURPRISE)]);
end

%%
% Assume that the the probability of drawing an a,b,c, or d are
% all equal to 0.25. Then the number of a's (or b's or c's or d's)
% drawing has a multnomial distribution with p(a) = p(b) = p(c) = p(d) =
% 0.25.

% Let's figure our what the correlation coefficient should be between
% any pair of spike triggered sums. We know it's going to be slightly
% negative because if you have an unusually large number of "a"s you are
% more likely to have a relativelysmall number of "b"s

% First theory
p1= .25; 
p2 = .25;
n = 100;  % number of spikes
nrepeats = 50;  % number of stixels

c = -n*p1*p2; % From Wikipedia
v1 = n*p1*(1-p1); % From Wikipedia
v2 = n*p2*(1-p2); % From Wikipedia
ro = c/sqrt(v1*v2); % <--- this is the correlation coefficient

% Here's a simulation that confirms the calculation above.
niter = 4000;
data =[];
for i = 1:niter
    a = unidrnd(4,n,nrepeats);
    b = [sum(a==1); sum(a==2)];
    % plot(b(1,:),b(2,:),'k.')
    r = corrcoef(b');
    data(i) = r(1,2);
end
figure; axes; hold on;
[counts,x] = hist(data);
bar(x,counts);
mean(data)

% Figuring out the sampling distribution of the correlation coefficient is
% too tricky for me to figure out how to do it. Instead, we'll pretend that
% the distribution of the numbers of "a"s and "b"s have a Gaussian
% distribution. Then we can use Fisher's Z transform to get the approximate
% sampling distribution.
x = linspace(-1,1,100);
z = 0.5*(log((1+x)./(1-x)) - log((1+ro)./(1-ro))).*sqrt(nrepeats-3)
y = diff(normcdf(z));
plot(x(1:end-1)+(x(2)-x(1))/2,(y./max(y))*max(counts),'m-','linewidth',2)
% Not too shabby for large n

%%
% Trying to get a confidence set for multinomial counts.
% In Patrick's experiment there are four types of stimuli: pLpM, pLmM, mLpM
% and mLmP. If there are 'n' spikes, we know that the number of stimuli in
% each of these categories has to sum to 'n'. This means we only have to
% think about the distribution of three of them ? the fourth doesn't provide
% any additional information. 

% Let's say we get a bunch of multinomial counts under the null hypothesis
% that the spikes have nothing to do with the stimulus. Where will these
% counts lie in the 3-D space?

alpha = 0.05; % significance level
n = 150;  % number of spikes
nrepeats = 100;  % number of stixels


% First I'm going to do a simulation
a = unidrnd(4,n,nrepeats);
b = [sum(a==1); sum(a==2); sum(a==3); sum(a==4)]; % simulating multinomial counts

figure; axes; hold on;
plot3(b(1,:),b(2,:),b(3,:),'.')

% I don't know how to make a 100*(1-alpha)% confidence set for a
% multinomial distribition. But I do know how to do it for a Gaussian
% distribution. Let's pretend the number of L-M's, L+M's, and -L-M's have a
% multivariate Gaussian distribution. To specify this distribution, I need
% to know the mean (3-D vector) and covariance matric (3 x 3 matrix).
% Luckily all the information I need is in Wikipedia.

% mu = [n/4 n/4 n/4]  (On average we expect n/4 spike-triggering sitmuli
% to be in any one of the categories.)

% the covariance matrix has n*p*(1-p) on the diagonal and -n.*p1*p2 on the
% off-diagnoal elements. These are just the variances and covariances of a
% multinomial ditribution (Wikipedia again).

% This code just makes the mean vector and covariance matrix
mu = n*[.25 .25 .25];
sigma = repmat(-n*.25^2,3,3);
sigma(1,1) = n*.25*.75; % This is n*p*(1-p)
sigma(2,2) = n*.25*.75;
sigma(3,3) = n*.25*.75;

% Sanity check stuff below - just making sure it works.
%[u,s,v] = svd(sigma);
%w = u*diag(1./sqrt(diag(s)));
%b_nomean = b-repmat(mean(b,2),1,size(b,2))
%tmp = w'*b_nomean([1 2 3],:)
%cov(tmp')

% OK, now using the fact that most draws from a Gaussian distribution fall
% close to the mean. How far from the mean has a chi-squared distribution.
crit = chi2inv(1-alpha,3);

% Doing significance test on each point (stixel)
ps = [];
for i = 1:nrepeats
    X = (b([1 2 3],i)-mean(b(1,:)))'*inv(sigma)*(b([1 2 3],i)-mean(b(1,:)));
    ps(i) = 1-chi2cdf(X,3);
end
L = ps<alpha;

sum(L)/nrepeats  % Here's the fraction of stixels that were found to be significant.
% It should be close to alpha.

% Doing some plotting just because it's fun to look at this stuff
% graphically.
[xx yy zz] = meshgrid(linspace(min(b(:)), max(b(:)), 10),...
    linspace(min(b(:)), max(b(:)), 10),...
    linspace(min(b(:)), max(b(:)), 10));

xformedxyz = [xx(:)-mean(b(1,:)) yy(:)-mean(b(2,:)) zz(:)-mean(b(3,:))];
X =[];
for i = 1:size(xformedxyz,1)
    X(i) = xformedxyz(i,:)*inv(sigma)*xformedxyz(i,:)';
end
surfstruct = isosurface(xx,yy,zz,reshape(X,size(xx)),crit);
h = patch(surfstruct);
set(h,'EdgeAlpha',.5,'FaceAlpha',.5,'FaceColor',[0 .5 0],'EdgeColor','none')
plot3(b(1,L),b(2,L),b(3,L),'r*')

% How about the sum of 'k' stixels?
n = 150;  % number of spikes
nrepeats = 100;  % number of total stixels
k = 4;  % number of selected stixels

% First I'm going to do a simulation
a = unidrnd(4,n,k,nrepeats);
b = [sum(a==1); sum(a==2); sum(a==3); sum(a==4)]; 
a = squeeze(sum(a,2))


%%
% finding points on an ellipse with equally spaced thetas
xmax = .09; ymax = .5;
thetas = linspace(0,2*pi,100);

% Wikipedia: http://en.wikipedia.org/wiki/Ellipse#Polar_form_relative_to_center
r = xmax*ymax./sqrt((ymax.*cos(thetas)).^2+(xmax.*sin(thetas)).^2);
polar(thetas,r,'.')

%%
% Finding the orientation and "spatial frequency" of one of Patrick's
% "tetris pieces" (from GridLMSubunit)


xys = [4 5; 4 6; 5 5; 5 6]; % a square
xys = [4 5; 4 6; 5 5; 5 6; 5 7]; % a square with an appendage
%xys = [4 5; 4 6; 4 7]; % a line
%xys = [1 1; 1 2; 2 2]; % an angle
%xys = [1 1; 1 2; 2 2;3 2]; % an "l"
%xys = [1 1; 1 2; 2 2;2 3]; % an "s"
%xys = [1 1; 2 2]; % diagonal
%xys = [4 4];  % a single stixel

if (size(xys,1) < 2)
    PC_1 = [1 0];  % single pixel, pick horizonatal
else
    coeff = pca(xys);
    PC_1 = coeff(:,1);
end
PC_end = [PC_1(2); -PC_1(1)]; 
theta = mod(atan2(PC_1(2),PC_1(1)),pi); % Preferred direction. 0 means horizonal

% The brute force way of getting the preferred spatial frequency:
% Getting all vertices
vertices = [];
canonicalvertices = [.5 .5; .5 -.5; -.5 -.5; -.5 .5];
for i = 1:size(xys,1)
    vertices = [vertices; repmat(xys(i,:),4,1)+canonicalvertices];
end
proj = vertices*PC_end;
halfperiod = range(proj);

% For computing the width of the projection of the tetris piece onto the
% axis that is orthogonal to the long axis:
% modtheta = mod(theta+pi/2,pi/2);
% if (modtheta <=pi/4)
%    squareprojwidth = cos(modtheta)+sin(modtheta).*tan(modtheta);
% else
%    squareprojwidth = sin(modtheta)+cos(modtheta).*cot(modtheta);
% end
% %squareprojwidth = cos(theta+pi/2)+sin(theta+pi/2);
% % projecting onto the vector orthogonal to the PC1
% proj = xys*PC_end;
% halfperiod = range(proj)+squareprojwidth % In stixels
%[theta halfperiod]
%sf = 1/(2*halfperiod); % cycles/stixel


% Plotting
figure; axes; hold on;
for i  = 1:size(xys,1)
    patch([-1 1 1 -1 -1]./2+xys(i,1),[-1 -1 1 1 -1]./2+xys(i,2),'red')
end
axis equal
[x,y] = pol2cart(theta+pi/2,halfperiod/2);
plot([-x x]+mean(xys(:,1)),[-y y]+mean(xys(:,2)),'g-','LineWidth',2);

% % Debugging
% angles = linspace(0,pi/2,100);
% % Zack's formula
% modangles = mod(angles,pi/2);
% squareprojwidth = zeros(size(angles));
% L = modangles<=pi/4;
% %squareprojwidth(L) = cos(modangles(L))+sin(modangles(L)).*tan(modangles(L));
% %squareprojwidth(~L) = sin(modangles(~L))+cos(modangles(~L)).*cot(modangles(~L));
% squareprojwidth = cos(modangles)+sin(modangles);
% plot(angles,squareprojwidth)


%%
% Plotting some GridLMSubunit data for a presentation
% Need to load a file first (e.g. LMSpikes_N052915001.mat).

% Evaluate this if using N052915001
LMSpikeMat = LMSpikes_N052915001;
LMSpikeMat = LMSpikeMat(1:150,:);
% Done with that

% Getting rid of high contrast LUM stimuli
L = softEq(LMSpikeMat(:,1), LMSpikeMat(:,2));
sqnorms = LMSpikeMat(:,1).^2+LMSpikeMat(:,1).^2;
LMSpikeMat(L & sqnorms > .15,:) = [];

%S = cov(LMSpikeMat(:,[1 2]));
%[u,s,v] = svd(S)
%wtmat = u*inv(sqrt(s));
%xy = LMSpikeMat(:,[1 2])*wtmat*v;
xy = LMSpikeMat(:,[1 2]);
%cov(xy) % sanity check
uniquestim = unique(LMSpikeMat(:,[1 2]),'rows');

maxmn = 0;
for i = 1:size(uniquestim,1)
    L = LMSpikeMat(:,1) == uniquestim(i,1) &  LMSpikeMat(:,2) == uniquestim(i,2);
    maxmn = max(maxmn, mean(LMSpikeMat(L,3)));
end

figure; axes; hold on;
for i = 1:size(uniquestim,1)
    L = LMSpikeMat(:,1) == uniquestim(i,1) &  LMSpikeMat(:,2) == uniquestim(i,2);
    mn = mean(LMSpikeMat(L,3));
    h = plot(xy(find(L,1),1),xy(find(L,1),2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([50*mn/maxmn 0]+4),'MarkerEdgeColor','white')
end
axis square;
set(gca,'Xlim',[-0.3 0.3],'Ylim',[-0.3 0.3]);

%%
% Plotting STA - all channels combined
M = [0.0608 0.1219 0.0175;
     0.0220 0.1266 0.0257;
     0.0019 0.0095 0.0976]; % From Dell4BitsCal(5)
 
GLMSPopAnalysis; % after running this line, click on the cell you want and then click on "Difference maps"
global DN
%whichframe = 5;
nstixperside = sqrt(size(DN.stats.pLpM.STA,1)/3);
STA = [];
% order of columns
%STAs = cat(3,DN.stats.pLpM.STA,DN.stats.mLmM.STA,DN.stats.pLmM.STA,DN.stats.mLpM.STA);
LMSTA = [];
LMSTA(:,:,1) = DN.stats.pLpM.STA(1:450,:)-.25;
LMSTA(:,:,2) = DN.stats.mLmM.STA(1:450,:)+.25;
LMSTA(:,:,3) = [DN.stats.pLmM.STA(1:225,:)-.25; DN.stats.pLmM.STA(226:450,:)+.25];
LMSTA(:,:,4) = [DN.stats.mLpM.STA(1:225,:)+.25; DN.stats.mLpM.STA(226:450,:)-.25];

weights = sum(sum(LMSTA.^2,1),3);
weights = weights-min(weights);
%STAs([451:675],:,:) = .25; % Setting S-cones to mean
%STA(:,1) = DN.stats.pLpM.STA(:,whichframe);
%STA(:,2) = DN.stats.mLmM.STA(:,whichframe);
%STA(:,3) = DN.stats.pLmM.STA(:,whichframe);
%STA(:,4) = DN.stats.mLpM.STA(:,whichframe);

greg = sum(LMSTA,3)*weights';
greg([451:675]) = 0; % adding zero s-cone components

lmsSTA = reshape(greg,nstixperside*nstixperside,3);
rgbSTA = lmsSTA*M; % For lms STA rgb = M'*[l m s]'
scalefactor = 1./(2*max(abs(rgbSTA(:))));
STAimg = scalefactor*reshape(rgbSTA,[nstixperside,nstixperside,3])+.5;
figure;
image(STAimg);
set(gca,'Xtick',[],'Ytick',[]);
axis square;

global GLMP
grid = [GLMP.subunit{1}.gridX{1}; GLMP.subunit{1}.gridY{1}];

%%
% Working on analyzing Patrick's psychometric data.
% First, preprocessing to isolate stimuli in each spoke
theta = atan2(gregmat(:,2),gregmat(:,1));
edges = linspace(-pi,pi,400);
n = histc(theta,edges);
if (n(end) > 0)
    n(end+1) = 0;
    edges(end+1) = edges(end) + edges(end)-edges(end-1)
end
bins = [edges(find(n>0))', edges(find(n>0)+1)']; % L and R edges of bins
spokeidx = nan(size(gregmat,1),1);
nstimperspoke = nan(size(bins,1),1);
thetas = nan(size(bins,1),1);
nspokes = size(bins,1);
for i = 1:nspokes
   L = theta>=bins(i,1) & theta<bins(i,2);
   spokeidx(L) = i; 
   nstimperspoke(i) = sum(L); % sanity check to make sure there are approx the same number of stim per spoke
   thetas(i) = bins(i,1);
end

% Now going spoke by spoke, fitting Weibull functions.
figure;
data = [];
for i = 1:nspokes
    L = spokeidx == i;
    r = sqrt(gregmat(L,1).^2+gregmat(L,2).^2)
    mn = gregmat(L,3);
    se = sqrt((mn.*(1-mn))./gregmat(L,4));
    subplot(ceil(sqrt(nspokes)),ceil(sqrt(nspokes)),i); hold on;
    errorbar(r,mn,se); 
    % Fitting 1-D weibull functions
    % "data" input argument should be [nCor nInc]? Yes.
    corinc = round(mn.*gregmat(L,4));
    corinc(:,2) = gregmat(L,4)-corinc(:,1);
    guesses = [mean(r) 2 max(mn)];
    [alpha, beta, gamma, success, modErrs] = weibullFit(r, corinc, 'mle', guesses);
    contrasts = linspace(min(r),max(r),100);
    plot(contrasts, gamma +(0.5-gamma).*exp(-((contrasts./alpha).^beta)));
    set(gca,'ylim',[.4 1]);
    data = [data; thetas(i) alpha beta gamma]
end

% figure;
% subplot(3,1,1); plot(data(:,1),data(:,2),'k-'); title('alpha');% alpha
% set(gca,'Xlim',[thetas(1) thetas(end)]);
% subplot(3,1,2); plot(data(:,1),data(:,3),'k-'); title('beta');% beta
% set(gca,'Xlim',[thetas(1) thetas(end)]);
% subplot(3,1,3); plot(data(:,1),data(:,4),'k-'); title('gamma');% gamma
% set(gca,'Xlim',[thetas(1) thetas(end)]);
% xlabel('theta (rad)');

% Now fitting alpha, beta, and gamma with ellipses
x = linspace(0,2*pi,100);
figure;
subplot(2,2,1); hold on;
polar(data([1:nspokes,1],1),data([1:nspokes,1],2)); % alpha
[tmpx,tmpy] = pol2cart(data(:,1), data(:,2));
tmpfn = @(params)ellipsefiterr(params,[tmpx tmpy],zeros(length(tmpx),1));
alpha_b = fmincon(tmpfn,[max(data(:,2)), min(data(:,2)),pi/4],[],[],[],[],[.01 .01 0],[10 10 2*pi])
polar(x+alpha_b(3),(alpha_b(1)*alpha_b(2))./sqrt((alpha_b(2)*cos(x)).^2+((alpha_b(1)*sin(x)).^2)))
title('alpha');

subplot(2,2,2); hold on;
polar(data([1:nspokes,1],1),data([1:nspokes,1],3)); % beta 
[tmpx,tmpy] = pol2cart(data(:,1), data(:,3));
tmpfn = @(params)ellipsefiterr(params,[tmpx tmpy],zeros(length(tmpx),1));
beta_b = fmincon(tmpfn,[max(data(:,2)), min(data(:,2)),pi/4],[],[],[],[],[.01 .01 0],[10 10 2*pi])
polar(x+beta_b(3),(beta_b(1)*beta_b(2))./sqrt((beta_b(2)*cos(x)).^2+((beta_b(1)*sin(x)).^2)))
title('beta');

subplot(2,2,3); hold on;
polar(data([1:nspokes,1],1),data([1:nspokes,1],4)); % gamma
[tmpx,tmpy] = pol2cart(data(:,1), data(:,3));
gamma_b = geomean(data(:,4));
polar(x,gamma_b*cos(x).^2+gamma_b*sin(x).^2)
title('gamma');
axis equal


% OK, now looking at all the fits using the Weibull parameters
% that vary ellipsoidally around the LM plane (except gamma which is fixed)
% Keep in mind that we have *not* fit the grand model to all of the data,
% we are fitting parameters of the individual spoke fits. This is just a
% diagnostic, not a final product.
figure;
for i = 1:nspokes
    L = spokeidx == i;
    r = sqrt(gregmat(L,1).^2+gregmat(L,2).^2);
    mn = gregmat(L,3);
    se = sqrt((mn.*(1-mn))./gregmat(L,4));
    subplot(ceil(sqrt(nspokes)),ceil(sqrt(nspokes)),i); hold on;
    errorbar(r,mn,se);
    
    alpha = (alpha_b(1)*alpha_b(2))./sqrt((alpha_b(2)*cos(thetas(i)-alpha_b(3))).^2+((alpha_b(1)*sin(thetas(i)-alpha_b(3))).^2))
    beta = (beta_b(1)*beta_b(2))./sqrt((beta_b(2)*cos(thetas(i)-beta_b(3))).^2+((beta_b(1)*sin(thetas(i)-beta_b(3))).^2))
    gamma = gamma_b;
    contrasts = linspace(min(r),max(r),100);
    plot(contrasts, gamma +(0.5-gamma).*exp(-((contrasts./alpha).^beta)));
    set(gca,'ylim',[.4 1]);
end

% Trying to do the whole fit at once
data = [];
for i = 1:nspokes
    L = spokeidx == i;
    mn = gregmat(L,3);
    cor = round(mn.*gregmat(L,4));
    inc = gregmat(L,4)-cor;
    data = [data; gregmat(L,1) gregmat(L,2) cor inc];
end
paramguess = [alpha_b(1) alpha_b(2) beta_b(1) beta_b(2) gamma alpha_b(3)]
tmpfn = @(params)WeibulSurfaceFiterr(params,data(:,[1 2]),data(:,[3 4]));
params = fminsearch(tmpfn,paramguess);

% Let's see how we did
figure;
for i = 1:nspokes
    L = spokeidx == i;
    r = sqrt(gregmat(L,1).^2+gregmat(L,2).^2);
    mn = gregmat(L,3);
    se = sqrt((mn.*(1-mn))./gregmat(L,4));
    subplot(ceil(sqrt(nspokes)),ceil(sqrt(nspokes)),i); hold on;
    errorbar(r,mn,se);
    
    alpha = (params(1)*params(2))./sqrt((params(2)*cos(thetas(i)-params(6))).^2+((params(1)*sin(thetas(i)-params(6))).^2))
    beta = (params(3)*params(4))./sqrt((params(4)*cos(thetas(i)-params(6))).^2+((params(3)*sin(thetas(i)-params(6))).^2))
    gamma = params(5);
    contrasts = linspace(min(r),max(r),100);
     plot(contrasts, gamma +(0.5-gamma).*exp(-((contrasts./alpha).^beta)));
    set(gca,'ylim',[.4 1]);
end

%%
% What does ridge regression look like wrt the interpretation of regression
% as  whitening the stimulus, calculating the response-weighted average
% and un-whitening the average by the "mechanism transform".

% ---------------
% Parameters you can manipulate
n = 100;
niter = 100;
noisesigma = 1;
mixmat = [1 .9; .9 1];
lambda = 1;
theta = pi/2;
trueb = [cos(theta) sin(theta)]';
% ---------------
paramestimates = zeros(niter,4);
for i = 1:niter
    X = normrnd(0,1,n,2)*mixmat;
    y = X*trueb+normrnd(0,noisesigma,n,1);
    b_est1 =(X'*X)\(X'*y);
%    b_est2 = inv(X'*X)*X'*y;
    % Now adding a ridge
    b_est2 =((X'*X)+eye(2)*lambda)\(X'*y);
    paramestimates(i,:) = [b_est1', b_est2'];
end
% Let's see how we did
figure; axes; hold on;
plot(paramestimates(:,1),paramestimates(:,2),'ks','MarkerFaceColor','black');
plot(paramestimates(:,3),paramestimates(:,4),'go','MarkerFaceColor','green');
%plot(trueb(1),trueb(2),'m*')
plot(0, 0,'m*');
plot([0 trueb(1)],[0 trueb(2)],'m-')
plot([0 mean(paramestimates(:,1))],[0 mean(paramestimates(:,2))],'k-')
plot([0 mean(paramestimates(:,3))],[0 mean(paramestimates(:,4))],'g-')
SSE1 = sum(sum((paramestimates(:,[1 2])-repmat(trueb',niter,1)).^2))
SSE2 = sum(sum((paramestimates(:,[3 4])-repmat(trueb',niter,1)).^2))
title(['SSE(OLS) = ',num2str(round(SSE1)),'     SSE(ridge) = ',num2str(round(SSE2))])

% Looking at the geometry
figure; subplot(2,2,1); hold on;
for i = 1:n
    h = plot(X(i,1),X(i,2),'ko');
    set(h,'MarkerSize',10*(y(i)-min(y))./(max(y)-min(y))+1)
end
whtmat1 = sqrtm(inv(X'*X));
whtmat2 = sqrtm(inv(X'*X+eye(2)*lambda));
whtmat2 = whtmat2*sqrt(det(whtmat1))/sqrt(det(whtmat2)); %equating variances of stimuli for display
% Remember that whtmat2 isn't really a whitening matrix
xformedX1 = X*whtmat1;
xformedX2 = X*whtmat2;
% sanity checks
cov(xformedX1);
cov(xformedX2);
det(whtmat1);
det(whtmat2);

subplot(2,2,2); hold on;
for i = 1:n
    h = plot(xformedX1(i,1),xformedX1(i,2),'ko');
    set(h,'MarkerSize',10*(y(i)-min(y))./(max(y)-min(y))+1)
    h = plot(xformedX2(i,1),xformedX2(i,2),'go');
    set(h,'MarkerSize',10*(y(i)-min(y))./(max(y)-min(y))+1)
end
% Finshing out the regression calculations just to make sure we're 
% doing everything correctly.
RWA1 = xformedX1'*y;
RWA2 = xformedX2'*y;
% bringing them back to the original space
% via the "mechanism transform". b_est1 and b_est2 are 
% still around from the iterated simulation, above.
b_est1_1 = RWA1'*whtmat1'
[b_est1_1;b_est1'] % worked!

b_est2_1 = RWA2'*whtmat2';
b_est2_1 = b_est2_1*norm(b_est2)./norm(b_est2_1);
[b_est2_1;b_est2']
mean(paramestimates)

%%
% Verifying that the STA/the variance of the stimulus dist (which is 
% identical in all directions because we're assuming radial symmetry)
% is the OLS estimate. Remember, variance has to be (sum of squares)/N, 
% not (sum of squares)/(N-1) for this to work.

sigma = 4;
nobs = 40;
x = sigma*[cos(linspace(0,2*pi,nobs+1));sin(linspace(0,2*pi,nobs+1))]';
x(end,:) =[];
% x, the stimulus distribution, is RS but does not have equal variance in
% all directions (because of the multiplication by sigma)

% Modeling a "half-squaring neuron"
weights = normrnd(0,1,2,1);
linresp= x*weights;
resp = max(0,linresp).^2

b_ols = regress(resp,x);
RWA = (resp'*x)/nobs;
RWS = (resp'*x);
b_rwa = RWA/var(x(:,1),1);

[b_ols'; b_rwa; weights']

%%
% OLS estimate of weighting vector magnitude is biased

niter = 10000;
b = [1 1 1];
nstim = 10;
sigma = .1;
data = zeros(niter, length(b));
for i = 1:niter
   x = normrnd(0,1,nstim, 3);
   y = x*b';
   y = y+normrnd(0,sigma,nstim,1);
   bhat = regress(y,x);
   data(i,:) = bhat';
end

norms = sqrt(sum(data.^2,2));
hist(norms);
realnorm = sqrt(b*b');
[h,p] = ttest(norms-realnorm)

%%
% Proving to myself that isoresponse contours of (stim*v1')^2+(stim*v2')=r
% are always elliptical, irrespective of v1 and v2 (they do not need to be
% orthogonal. (See GrantBrainStorming Section 4.10 for a continuation of
% this exercise)

[x,y] = meshgrid([-10:.1:10],[-10:.1:10]);
v1 = [1./sqrt(2); -1./sqrt(2)];
v2 = [1; 1];

fv = [x(:) y(:)]*[v1 v2];
fv = sum(abs(fv).^2,2);
fv = reshape(fv,size(x));
figure; 
contour(x,y,fv);

% Finding a pair of (v1,v2) for which the function (stim*v1')^2+(stim*v2')=r
% is identical.
% (stim*v1')^2+(stim*v2')=1
% ([L M]*[a b]')^2+([L M]*[c d]')^2=1
% (a^2+c^2)L+2(ab+cd)LM+(b^2+d^2)M=1
% Let
% A = a^2+c^2
% B = ab+cd
% C = b^2+d^2

myfun = @(x) [x(1)^2 + x(3)^2;x(1)*x(2)+x(3)*x(4); x(2)^2+x(4)^2]
% starting with this ellipse
a = v1(1); b = v1(2);
c = v2(1); d = v2(2);
ABC = myfun([a b c d]);

myerrfun = @(x)sum((myfun(x)-ABC).^2);
% testing
myerrfun(normrnd(0,1,4,1))
myerrfun([a b c d])

niter = 100;
data = [];
for i = 1:niter
    newabcd = fminsearch(myerrfun,normrnd(0,1,4,1));
    fv0 = [x(:) y(:)]*[v1 v2];
    fv0 = sum(abs(fv).^2,2);
    fv1 = [x(:) y(:)]*reshape(newabcd,2,2);
    fv1 = sum(abs(fv).^2,2);
    if ~all(fv1 == fv0)
        keyboard
    end
    data = [data, newabcd];
end


