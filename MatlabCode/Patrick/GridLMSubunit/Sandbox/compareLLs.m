global GLMSPopData

desinfo = GLMSPopData(1,:);
LLx10 = cat(1,GLMSPopData{2:end,strcmp(desinfo,'LLx10')});
LLx2 = cat(1,GLMSPopData{2:end,strcmp(desinfo,'LLx2')});
whichnd = GLMSPopData(2:end,strcmp(desinfo,'nD'))
oneDL = strcmp(whichnd,'1D')
twoDL = strcmp(whichnd,'2D')
surfparams = cat(1,GLMSPopData{2:end,strcmp(desinfo,'Surface Parameters')})
surfparams = surfparams(2:2:end);
oneDparams = surfparams(oneDL);
twoDparams = cat(1,surfparams{twoDL});


LLx101D = cat(1,LLx10.oneD);
LLx102D = cat(1,LLx10.twoD);
LLx21D = cat(1,LLx2.oneD);
LLx22D = cat(1,LLx2.twoD);


figure; clf; hold on;
plot(LLx101D,LLx21D,'bo')
plot(LLx102D,LLx22D,'ro')
xlabel('2D at n=10')
ylabel('2D at n=2')
plot([0, 18000],[0,18000],'k--')