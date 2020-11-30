function dispimage(~,~,x,y,basis,bkgnd)
im = x*(basis{1}-bkgnd)+y*(basis{2}-bkgnd);
im = 0.5 + 0.5*im/(max(abs(im(:)))+0.01);
figure;image(im); axis square; set(gca,'XTick',[],'YTick',[]);
end