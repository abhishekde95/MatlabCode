% Conway and Livingstone, 2006 stim design and simulated DO response
% Author - Gregory D. Horwitz, 1/20

close all; clearvars;
x = -6:6;
RF = normpdf(x,0,1)'*normpdf(x,0,1);
imagesc(RF)
 
% Simulated experiment in which a single pixel increment and decrement 
% of light appear on every frame. 
% Also assuming the the response of the cell is purely linear.
STS = zeros(size(RF));
niter = 10000;
for i = 1:niter
    im = zeros(size(RF));
    inc_pos = unidrnd(size(RF)); % position of increment
    dec_pos = unidrnd(size(RF)); % position of decrement
    % setting up image
    im(inc_pos(2),inc_pos(1)) =  im(inc_pos(2),inc_pos(1)) + 1;
    im(dec_pos(2),dec_pos(1)) =  im(dec_pos(2),dec_pos(1)) - 1;
    STS = STS+im*(im(:)'*RF(:));
end

STSim = 0.5 + (0.5*(STS - mean(STS(:)))/(max(abs(STS(:)))+eps));
figure; 
subplot(2,2,1);
imagesc(STSim); axis square; axis xy; set(gca,'Tickdir','out');
subplot(2,2,2); hold on;
plot(STS(:)); 
plot([1 numel(RF)],[0 0],'-'); axis square; set(gca,'Tickdir','out');
 
% Trying again with a spiking nonlinearity
STS = zeros(size(RF));
for i = 1:niter
    im = zeros(size(RF));
    inc_pos = unidrnd(size(RF)); % position of increment
    dec_pos = unidrnd(size(RF)); % position of decrement
    % setting up image
    im(inc_pos(2),inc_pos(1)) =  im(inc_pos(2),inc_pos(1)) + 1;
    im(dec_pos(2),dec_pos(1)) =  im(dec_pos(2),dec_pos(1)) - 1;
    STS = STS+im*max(0,im(:)'*RF(:)); % Here's the nonlinearity
end
STSim = 0.5 + (0.5*(STS - mean(STS(:)))/(max(abs(STS(:)))+eps));
subplot(2,2,3);
imagesc(STSim);  axis square; axis xy; set(gca,'Tickdir','out');
subplot(2,2,4); hold on;
plot(STS(:));
plot([1 numel(RF)],[0 0],'-'); axis square;  set(gca,'Tickdir','out');