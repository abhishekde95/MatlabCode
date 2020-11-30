function [sigma, mu, A, yfit, xx, yyfit] = fitGaussian(x,y)

p = polyfit(x, log(y), 2);
sigma = sqrt(-1/(2*p(1)));
mu = p(2)*sigma^2;
A = exp(p(3)+mu^2/(2*sigma^2));
yfit = A * exp( -(x-mu).^2 / (2*sigma^2) );

%finer fit
xx = min(x): min(x)/100 :max(x); 
yyfit = A * exp( -(xx-mu).^2 / (2*sigma^2) );
end