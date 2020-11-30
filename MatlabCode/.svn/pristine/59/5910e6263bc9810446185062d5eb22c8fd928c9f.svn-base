% For playing with RWA

xrange = 20;
x = linspace(-xrange,xrange,20);
nsp = ComputeNakaRushtonJPW([100 mean(abs(x)) 3 0],x,'symmetric');
%stanorm = mean(nsp./sum(nsp).*x)./var(x,1)
staraw = mean(nsp.*x)./var(x,1)
b = regress(nsp',x')

figure(44); clf; hold on; box on;
% for n = 1:5
%     plot(x,poissrnd(nsp),'ko')
% end
plot(x,nsp,'r*')
plot(x,nsp,'r-')
plot(x,x*b,'b-')
%plot(x,x*stanorm,'b--')
plot(x,x*staraw,'g--')

% %% Format and save figure
% set(gcf,'PaperPositionMode','auto')
% if ismac
%     print('/Users/jpatrickweller/Dropbox/VisionResearchPaper/Fig1','-depsc')
% else ispc
%     print('C:\Users\jpweller\Dropbox\VisionResearchPaper\Fig1','-depsc');
% end
% disp('Fig 1 done.')