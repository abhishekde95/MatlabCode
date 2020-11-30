function [ori_tuning, sf_tuning] = resp_to_sinusoids(out,lambda,theta,mode)
% Extracting sf and orituning

N = 21;
[X, Y] = meshgrid(1:1:N,1:1:N);
med = ceil(median(1:N));
X = X-med; Y = Y-med; % So negative numbers mean down

resp= zeros(numel(lambda),numel(theta));
if mode == 1 % Crescent
    X1 = X-out.muxc; Y1 = Y+out.muyc; % For center, So negative numbers mean down
    X2 = X-out.muxs; Y2 = Y+out.muys; % For surround, So negative numbers mean down
    im = out.gaincenter*sqrt(1/(2*pi*out.centersigma^2))*exp(-(X1.^2+Y1.^2)/out.centersigma^2) - out.gainsurround*sqrt(1/(2*pi*out.surroundsigma^2))*exp(-(X2.^2+Y2.^2)/out.surroundsigma^2);
elseif mode == 2 % Gabor 
    xprime = X.*cos(-out.theta)+Y.*sin(-out.theta);
    yprime = -X.*sin(-out.theta)+Y.*cos(-out.theta);
    im = out.amplitude*exp(-(xprime.^2+out.gamma.^2.*yprime.^2)./(2.*out.sigma.^2)).*cos((2.*pi.*yprime./out.lambda)-out.phi);
elseif mode == 3 % DOG
    im = out.gaincenter*sqrt(1/(2*pi*out.centersigma^2))*exp(-(X.^2+Y.^2)/out.centersigma^2) - out.gainsurround*sqrt(1/(2*pi*out.surroundsigma^2))*exp(-(X.^2+Y.^2)/out.surroundsigma^2);
end

count = 1;
for ii = 1:numel(lambda)
    for jj = 1:numel(theta)
        xprime = X.*cos(-theta(jj))+Y.*sin(-theta(jj));
        yprime = -X.*sin(-theta(jj))+Y.*cos(-theta(jj));
        sinusoid = cos((2.*pi.*yprime./lambda(ii)));
        resp(ii,jj) = abs(real(sum(sum(fftshift(fft2(sinusoid)).*fftshift(fft2(im))))));
%         figure(1),subplot(numel(lambda),numel(theta),count),imagesc(sinusoid), set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray');
%         figure(2),subplot(numel(lambda),numel(theta),count),imagesc(abs(real(fftshift(fft2(sinusoid))))), set(gca,'XTick',[],'YTick',[]); axis square; colormap('gray');
        count = count + 1;
    end
end
ori_tuning = sum(resp,1); 
ori_tuning = ori_tuning/max(ori_tuning);
sf_tuning = sum(resp,2); 
sf_tuning = sf_tuning/max(sf_tuning);
end

