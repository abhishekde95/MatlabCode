% 2D slices

% clear;
% clc;

% load 5thCellAlldata.mat
% load 2ndCellAlldata.mat
% load LMcellAlldata.mat

%%
lms = unique(support_x, 'rows');
lengAxis = 100;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

[ll, mm] = meshgrid(xaxis, yaxis);

% F = TriScatteredInterp(lms(:,1), lms(:,2), spk);
% F = TriScatteredInterp(lms(:,1), lms(:,2), nspikes);
% F = TriScatteredInterp(lms(:,1), lms(:,2), log(1+exp(fmap)));
% F = TriScatteredInterp(lms(:,1), lms(:,2), exp(fmapFinal));
F = TriScatteredInterp(lms(:,1), lms(:,2), g(predictiveMean));
qz = F(ll, mm); 
figure(1); 
% subplot(221); plot(lms(:,1), 'b.'); hold on; plot(lms(:,2), 'r.');
subplot(122); contour(ll, mm, qz);  title('Mean with one lengthscale'); ylabel('m-cone'); xlabel('l-cone');
% hold on; scatter(lms(:,1), lms(:,2), 'o');

%%
avg_r = (spk(1:3:end)+spk(2:3:end)+spk(3:3:end))./3; 
F2 = TriScatteredInterp(lms(:,1), lms(:,2), avg_r);
% F2 = TriScatteredInterp(lms(:,1), lms(:,2), log(exp(datastruct.finit)+1));
qz2 = F2(ll, mm); 
subplot(122); surf(ll, mm, qz2); title('avg r'); ylabel('m-cone'); xlabel('l-cone');


%%

clear;
clc;

% load 5thCellAlldata.mat
% load 4thCellAlldata.mat
load 1stCellAlldata.mat
% load 2ndCellAlldata.mat
lms = (support_x);
lengAxis = 20;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
numbBinz = length(zaxis)


% 2D slice

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:) zeros(length(ll(:)), 1)];
datastruct.support = llmm;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
figure(1);
subplot(5, 3, 1);
% subplot(221); 
contour(ll, mm, g(Mean_LM)); title('LM'); ylabel('cell 1');

[ll, ss] = meshgrid(xaxis, zaxis);
llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
datastruct.support = llss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 2); contour(ll, ss, g(Mean_LS)); axis xy; title('LS')

[mm, ss] = meshgrid(yaxis, zaxis);
mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
datastruct.support = mmss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 3); contour(mm, ss, g(Mean_MS)); axis xy; title('MS')

%%

clear;
clc;

% load 5thCellAlldata.mat
% load 4thCellAlldata.mat
% load 1stCellAlldata.mat
load 2ndCellAlldata.mat

lms = (support_x);
lengAxis = 20;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
numbBinz = length(zaxis)


% 2D slice

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:) zeros(length(ll(:)), 1)];
datastruct.support = llmm;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
figure(1);
subplot(5, 3, 4);
% subplot(221); 
contour(ll, mm, g(Mean_LM)); title('LM'); ylabel('cell 2');

[ll, ss] = meshgrid(xaxis, zaxis);
llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
datastruct.support = llss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 5); contour(ll, ss, g(Mean_LS)); axis xy; title('LS'); 

[mm, ss] = meshgrid(yaxis, zaxis);
mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
datastruct.support = mmss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 6); contour(mm, ss, g(Mean_MS)); axis xy; title('MS');

%%
clear;
clc;

% load 5thCellAlldata.mat
% load 4thCellAlldata.mat
load 3rdCellAlldata.mat
% load 2ndCellAlldata.mat
lms = (support_x);
lengAxis = 20;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
numbBinz = length(zaxis)


% 2D slice

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:) zeros(length(ll(:)), 1)];
datastruct.support = llmm;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
figure(1);
subplot(5, 3, 7);
% subplot(221); 
contour(ll, mm, g(Mean_LM)); title('LM'); ylabel('cell 3');

[ll, ss] = meshgrid(xaxis, zaxis);
llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
datastruct.support = llss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 8); contour(ll, ss, g(Mean_LS)); axis xy; title('LS');

[mm, ss] = meshgrid(yaxis, zaxis);
mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
datastruct.support = mmss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 9); contour(mm, ss, g(Mean_MS)); axis xy; title('MS');

%%

clear;
clc;

% load 5thCellAlldata.mat
load 4thCellAlldata.mat
% load 1stCellAlldata.mat

lms = (support_x);
lengAxis = 20;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
numbBinz = length(zaxis)


% 2D slice

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:) zeros(length(ll(:)), 1)];
datastruct.support = llmm;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
figure(1);
subplot(5, 3, 10);
% subplot(221); 
contour(ll, mm, g(Mean_LM)); title('LM'); ylabel('cell 4');

[ll, ss] = meshgrid(xaxis, zaxis);
llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
datastruct.support = llss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 11); contour(ll, ss, g(Mean_LS)); axis xy; title('LS');

[mm, ss] = meshgrid(yaxis, zaxis);
mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
datastruct.support = mmss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 12); contour(mm, ss, g(Mean_MS)); axis xy; title('MS');

%%

clear;
clc;

load 5thCellAlldata.mat

lms = (support_x);
lengAxis = 20;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
numbBinz = length(zaxis)


% 2D slice

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:) zeros(length(ll(:)), 1)];
datastruct.support = llmm;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
figure(1);
subplot(5, 3, 13);
% subplot(221); 
contour(ll, mm, g(Mean_LM)); title('LM'); ylabel('cell 5');

[ll, ss] = meshgrid(xaxis, zaxis);
llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
datastruct.support = llss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 14); contour(ll, ss, g(Mean_LS)); axis xy; title('LS');

[mm, ss] = meshgrid(yaxis, zaxis);
mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
datastruct.support = mmss;
predictiveMean = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(5, 3, 15); contour(mm, ss, g(Mean_MS)); axis xy; title('MS');

%%
% load 1stCellAlldata.mat
load 4thCellAlldata.mat

lms = (support_x);
lengAxis = 20;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
numbBinz = length(zaxis)

% 2D slice

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:) zeros(length(ll(:)), 1)];
datastruct.support = llmm;
[predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
figure(2);
subplot(333);
% subplot(221); 
contour(ll, mm, g(Mean_LM)); title('total')

[ll, ss] = meshgrid(xaxis, zaxis);
llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
datastruct.support = llss;
[predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(336); contour(ll, ss, g(Mean_LS)); axis xy; 

[mm, ss] = meshgrid(yaxis, zaxis);
mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
datastruct.support = mmss;
[predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
subplot(339); contour(mm, ss, g(Mean_MS)); axis xy; 


%% mse 

% load 4thCellAlldata.mat
% load 2ndCellAlldata.mat
load LMcellAlldata.mat

truelambda = exp(fmapFinal);

% load f_2nd_infomax_150;
% load f_2nd_rand_150;
% a = sum((truelambda- exp(f_2nd_infomax_150)).^2)/length(support_x)
% b = sum((truelambda- exp(f_2nd_rand_150)).^2)./length(support_x)

% load f_4th_infomax_150;
% load f_4th_rand_150;
% load f_ml_infomax_150;
% load f_ml_rand_150;
load f_ml_infomax_300;
load f_ml_rand_300;
a = sum((truelambda- exp(f_ml_infomax_300)).^2)/length(support_x)
b = sum((truelambda- exp(f_ml_rand_300)).^2)./length(support_x)

%%
% clear;
% clc;

% load 5thCellMaxInfo150.mat
% load 2ndCellMaxInfo150.mat
load mlCellrand300.mat;


F = TriScatteredInterp(lms(:,1), lms(:,2), exp(f_ml_rand_300));
qz = F(ll, mm); 
figure(1); subplot(224); contour(ll, mm, qz); axis image; title(['rand(300), mse:' num2str(b)])


%%

lms = (support_x);
lengAxis = 50;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

% zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
% numbBinz = length(zaxis)

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:)];
datastruct.support = llmm;
[predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
% figure(233);
subplot(132); contour(ll, mm, g(Mean_LM)); axis image;  
title(['rand(200), mse:' num2str(b)])

% [ll, ss] = meshgrid(xaxis, zaxis);
% llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
% datastruct.support = llss;
% [predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
% Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
% subplot(335); contour(ll, ss, g(Mean_LS));axis xy; 
% 
% [mm, ss] = meshgrid(yaxis, zaxis);
% mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
% datastruct.support = mmss;
% [predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
% Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
% subplot(338); contour(mm, ss, g(Mean_MS)); axis xy; 

%%
% clear;
% clc;

% load 2ndCellrand150.mat
% load 5thCellrand150.mat;
load mlCellMaxInfo300.mat;


F = TriScatteredInterp(lms(:,1), lms(:,2), exp(f_ml_infomax_300));
qz = F(ll, mm); 
figure(1); subplot(222); contour(ll, mm, qz); axis image; title(['info(300), mse:' num2str(a)]); 


%%

lms = (support_x);
lengAxis = 50;
xaxis = linspace(min(lms(:,1)), max(lms(:,1)), lengAxis); % for cell0
numbBinx = length(xaxis)

yaxis = linspace(min(lms(:,2)), max(lms(:,2)), lengAxis); % for cell0
numbBiny = length(yaxis)

% zaxis = linspace(min(lms(:,3)), max(lms(:,3)), lengAxis);
% numbBinz = length(zaxis)

[ll, mm] = meshgrid(xaxis, yaxis);
llmm = [ll(:) mm(:)];
datastruct.support = llmm;
[predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
Mean_LM = reshape(predictiveMean, numbBinx, numbBinx);
% figure(233);
subplot(133); contour(ll, mm, g(Mean_LM)); axis image; 
title(['info(200), mse:' num2str(a)]); 

% [ll, ss] = meshgrid(xaxis, zaxis);
% llss = [ll(:) zeros(length(ll(:)), 1) ss(:)];
% datastruct.support = llss;
% [predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
% Mean_LS = reshape(predictiveMean, numbBinx, numbBinx);
% subplot(334); contour(ll, ss, g(Mean_LS)); ylabel('l-s');axis xy; 
% 
% [mm, ss] = meshgrid(yaxis, zaxis);
% mmss = [ zeros(length(mm(:)), 1) mm(:) ss(:)];
% datastruct.support = mmss;
% [predictiveMean, predictiveVar, xNext, idxNext] = makePrediction(prs, datastruct, aFinal, WFinal, sqrtLFinal);
% Mean_MS = reshape(predictiveMean, numbBinx, numbBinx);
% subplot(337); contour(mm, ss, g(Mean_MS)); ylabel('m-s');axis xy; 

