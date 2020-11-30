% Simulating some experimental outcomes for the ScottMurray project.

npts = 400;
annuluswidth = .3;
[x,y] = meshgrid (linspace(-2,2,npts), linspace(-2,2,npts));
squareddist = x.^2+y.^2;
annulus = single(squareddist < 1 & squareddist > 1-annuluswidth);

RFsize = [10 50 100 200];
for i = 1:length(RFsize)
    x = normpdf(linspace(-2,2,RFsize(i)),0,1);
    RF = x'*x;
    
    blurredim = conv2(annulus, RF, 'same');
    subplot(2,length(RFsize),i);
    imagesc(blurredim);
    axis image;
    colormap(gray);
    set(gca,'YTick',[],'XTick',[]);
    subplot(2, length(RFsize),i+length(RFsize));
    plot(blurredim(round(npts/2),:));
    axis square;
    set(gca,'YTick',[]);
    peak = find(blurredim(round(npts/2),[1:round(npts/2)]) == max(blurredim(round(npts/2),[1:round(npts/2)])),1);
    title(['Peak is at x = ',num2str(peak)]);
end