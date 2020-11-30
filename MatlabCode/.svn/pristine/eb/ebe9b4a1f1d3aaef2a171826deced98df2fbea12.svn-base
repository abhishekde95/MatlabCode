% A script for implememting robust regression using tukey bisquare error
% Author - Abhishek De, 12/19
close all; clearvars;

e = linspace(-10,10,1001);

y_tb = []; w_tb = []; k_tb = 4.685; % Tukey Bisquare
y_hu = []; w_hu = []; k_hu = 1.345; % Huber 
for ii = 1:numel(e)
    if abs(e(ii))<k_tb
        y_tb = [y_tb; (k_tb^2/6)*(1-(1-(e(ii)/k_tb).^2).^3)];
        w_tb = [w_tb; (1-(e(ii)/k_tb).^2).^2];
    elseif abs(e(ii))>=k_tb
        y_tb = [y_tb; k_tb.^2/6];
        w_tb = [w_tb; 0];
    end
    
    if abs(e(ii))<k_hu
        y_hu = [y_hu; (e(ii).^2)/2];
        w_hu = [w_hu; 1];
    elseif abs(e(ii))>=k_hu
        y_hu = [y_hu; k_hu*abs(e(ii))-0.5*k_hu.^2];
        w_hu = [w_hu; k_hu/abs(e(ii))];
    end
end

figure(1); subplot(121); plot(e,y_tb,'k'); hold on; plot(e,w_tb,'r'); plot(e,y_tb.*w_tb,'g'); axis square; set(gca,'Tickdir','out'); xlabel('error'); title('Tukey-Bisquare');
subplot(122); plot(e,y_hu,'k'); hold on; plot(e,w_hu,'r'); plot(e,y_hu.*w_hu,'g'); axis square; set(gca,'Tickdir','out'); xlabel('error'); title('Huber');