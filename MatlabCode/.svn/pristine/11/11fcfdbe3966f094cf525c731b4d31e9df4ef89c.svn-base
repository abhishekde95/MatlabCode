% Analysis of Greg's ISORESPONSE data (Horwitz and Haas; 2012)
% Author - Abhishek De
close all; clearvars;
load GregISORESPONSE_data.mat
load eigenvalues.mat
load eigenvectors.mat
load p_ftest.mat % Comparing the linear and quadratic model 
plot_counter = 1;
% data = [WHICHSURF modulationratio prefSF coneweights]
% WHICHSURF=1 Hyperboloid of 1 sheet
% WHICHSURF=2 Hyperboloid of 2 sheet
% WHICHSURF=3 Ellipsoid

% Plotting the cone weights of cells described by planes\
SIMPLE_ID = data(:,2)>1;
PLANES_ID = p_ftest>=0.01;
conewts_planesS = data(PLANES_ID & SIMPLE_ID,4:6);
conewts_planesS = conewts_planesS.*repmat(sign(conewts_planesS(:,2)),[1 3]);
conewts_planesC = data(PLANES_ID & ~SIMPLE_ID,4:6);
conewts_planesC = conewts_planesC.*repmat(sign(conewts_planesC(:,2)),[1 3]);

%PLotting cone weights for cells best described by planes
figure(plot_counter); subplot(311); plot(conewts_planesS(:,1),conewts_planesS(:,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_planesC(:,1),conewts_planesC(:,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]);
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k'); 
xlabel('L'), ylabel('M'); legend('S','C'); title('Planes');

% Calculating the principal directions of the cells better fit by Hyperboloid of 1 or 2 sheet

for ii=1:size(data,1)
    if data(ii,1)==1 & p_ftest(ii) <0.01
        wts = data(ii,[4:6])'; wts = wts*sign(wts(2));
%         wts = eigenvectors{ii}(:,2);
%         wts = wts./repmat(sum(abs(wts),1),[3 1]);
%         wts = wts.*repmat(sign(wts(2,:)),[3 1]);
        if data(ii,2)> 1
            figure(plot_counter); subplot(312); plot(wts(1,:),wts(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        else
           figure(plot_counter); subplot(312); plot(wts(1,:),wts(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
        end
    elseif data(ii,1) == 2 & p_ftest(ii)<0.01
        wts = data(ii,[4:6])'; wts = wts*sign(wts(2));
%         wts = eigenvectors{ii}(:,3);
%         wts = wts./repmat(sum(abs(wts),1),[3 1]);
%         wts = wts.*repmat(sign(wts(2,:)),[3 1]);
        if data(ii,2)> 1
            figure(plot_counter); subplot(313); plot(wts(1,:),wts(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
        else
           figure(plot_counter); subplot(313); plot(wts(1,:),wts(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0]); hold on;
        end
    elseif data(ii,1) == 3
       
    end
end

figure(plot_counter); subplot(312); axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k'); 
xlabel('L'), ylabel('M'); title('Hyperboloiod of 1 sheet');
figure(plot_counter); subplot(313); axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k'); 
xlabel('L'), ylabel('M'); title('Hyperboloiod of 2 sheets');


