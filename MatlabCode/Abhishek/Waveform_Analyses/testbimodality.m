function [p,occurences] = testbimodality(h,x,xmin,xmax,nboot)

vec = xmin:h:xmax;
L = numel(vec);
N = numel(x);
d = zeros(L,1);
for ii = 1:L
    d(ii) = (1/(N*h))*sum(exp((-((vec(ii)-x)/h).^2/2)*sqrt(2)*pi));
end

% figure(1); plot(vec,d,'k','Linewidth',2); axis square; set(gca,'Tickdir','out');

nmodes = zeros(nboot,1);
for jj = 1:nboot
    d_new = zeros(L,1);
    x_new = randsample(xmin:h:xmax,N,true,d);
    for ii = 1:L
        d_new(ii) = (1/(N*h))*sum(exp((-((vec(ii)-x_new)/h).^2/2)*sqrt(2)*pi));
    end
    nmodes(jj) = sum(islocalmax(d_new));
end
p = sum(nmodes<=1)/nboot;

occurences = [sum(nmodes==1) sum(nmodes>1)];

end

