% Author - Abhishek De 
% Script for generating stimulus in cone space
function gen_color__stim()
clearvars;
close all;
load fundamentals.mat % CONE FUNDAMENTALS: L,M,S
load mon_spd.mat % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
 
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S 
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = mon_spd*fundamentals; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
invM = inv(M);
SetupFig(invM);

end

function SetupFig(in)
gcf = figure(1);
set(gcf,'DefaultAxesUnits','pixels')
set(gcf,'position',[10 10 1000 800]);
set(gcf,'ButtonDownFcn','drawnow');  % in case user drags figure
clf;
a = get(gcf,'UserData');
a.param.invM = in;

% Stimulus 1 features 
a.uicontrols.stim_label1 = uicontrol('style','text','String','STIM 1','Position',[100 330 60 20]);
a.uicontrols.textL1 = uicontrol('style','text','String','L:0.1','Position',[100 50 60 20]);
a.uicontrols.L1 = uicontrol('style','slider','Min',-1,'Max',1,'Value',0.1,'Position',[100 20 60 20],'Callback',@show_stim1);
a.uicontrols.textM1 = uicontrol('style','text','String','M:-0.1','Position',[170 50 60 20]);
a.uicontrols.M1 = uicontrol('style','slider','Min',-1,'Max',1,'Value',-0.1,'Position',[170 20 60 20],'Callback',@show_stim1);
a.uicontrols.textS1 = uicontrol('style','text','String','S:-0.6','Position',[240 50 60 20]);
a.uicontrols.S1 = uicontrol('style','slider','Min',-10,'Max',10,'Value',-0.6,'Position',[240 20 60 20],'Callback',@show_stim1);
a.uicontrols.textsf1 = uicontrol('style','text','String','sf:1','Position',[175 360 60 20]);
a.uicontrols.sf1 = uicontrol('style','slider','Min',0,'Max',10,'Value',1,'Position',[175 330 60 20],'Callback',@show_stim1);
a.uicontrols.textphase1 = uicontrol('style','text','String','phase:0','Position',[250 360 100 20]);
a.uicontrols.phase1 = uicontrol('style','slider','Min',0,'Max',pi,'Value',0,'Position',[250 330 100 20],'Callback',@show_stim1);
a.axeshandles.stim1 = axes('position',[100 100 200 200]);
set(a.axeshandles.stim1,'XTick',[],'YTick',[],'Box','on');
axis image;
a.stim1 = [];

% Stimulus 2 features 
a.uicontrols.stim_label2 = uicontrol('style','text','String','STIM 2','Position',[100 710 60 20]);
a.uicontrols.textL2 = uicontrol('style','text','String','L:0.1','Position',[100 430 60 20]);
a.uicontrols.L2 = uicontrol('style','slider','Min',-1,'Max',1,'Value',0.1,'Position',[100 400 60 20],'Callback',@show_stim2);
a.uicontrols.textM2 = uicontrol('style','text','String','M:-0.1','Position',[170 430 60 20]);
a.uicontrols.M2 = uicontrol('style','slider','Min',-1,'Max',1,'Value',-0.1,'Position',[170 400 60 20],'Callback',@show_stim2);
a.uicontrols.textS2 = uicontrol('style','text','String','S:0.6','Position',[240 430 60 20]);
a.uicontrols.S2 = uicontrol('style','slider','Min',-10,'Max',10,'Value',0.6,'Position',[240 400 60 20],'Callback',@show_stim2);
a.uicontrols.textsf2 = uicontrol('style','text','String','sf:1','Position',[175 740 60 20]);
a.uicontrols.sf2 = uicontrol('style','slider','Min',0,'Max',10,'Value',1,'Position',[175 710 60 20],'Callback',@show_stim2);
a.uicontrols.textphase2 = uicontrol('style','text','String','phase:0','Position',[250 740 100 20]);
a.uicontrols.phase2 = uicontrol('style','slider','Min',0,'Max',pi,'Value',0,'Position',[250 710 100 20],'Callback',@show_stim2);
a.axeshandles.stim2 = axes('position',[100 480 200 200]);
set(a.axeshandles.stim2,'XTick',[],'YTick',[],'Box','on');
axis image;
a.stim2 = [];

% Template features
theta = 0; sigma = 0.4; ggamma = 1;
[x,y] = meshgrid(linspace(-1,1,51), linspace(-1,1,51));
X = x*cos(-theta) + y*sin(-theta);
Y =-x*sin(-theta) + y*cos(-theta);
a.template.expterm = exp(-(X.^2 + ggamma^2*Y.^2)/2/sigma^2);
a.template.costerm = cos(2*pi*Y*1.0);
a.template.X = X;
a.template.Y = Y;

% S-OFF receptive field 
a.uicontrols.soff_label = uicontrol('style','text','String','S OFF','Position',[400 330 60 20]);
a.uicontrols.textsigma_soff = uicontrol('style','text','String','sigma:0.25','Position',[500 360 100 20]);
a.uicontrols.sigma_soff = uicontrol('style','slider','Min',0,'Max',2.0,'Value',0.25,'Position',[500 330 100 20],'Callback',@show_soffrf);
a.axeshandles.SOFF = axes('position',[400 100 200 200]);
set(a.axeshandles.SOFF,'XTick',[],'YTick',[],'Box','on');
axis image;
a.SOFF_rgb = [];

% S-ON receptive field 
a.uicontrols.son_label = uicontrol('style','text','String','S ON','Position',[700 330 60 20]);
a.uicontrols.textsigma_son = uicontrol('style','text','String','sigma:0.25','Position',[800 360 100 20]);
a.uicontrols.sigma_son = uicontrol('style','slider','Min',0,'Max',2.0,'Value',0.25,'Position',[800 330 100 20],'Callback',@show_sonrf);
a.axeshandles.SON = axes('position',[700 100 200 200]);
set(a.axeshandles.SON,'XTick',[],'YTick',[],'Box','on');
axis image;
a.SON_rgb = [];

% UICONTROL for jointly varying the size of receptive 
a.uicontrols.text_factor = uicontrol('style','text','String','factor:1.0','Position',[500 50 100 20]);
a.uicontrols.factor = uicontrol('style','slider','Min',0.1,'Max',2,'Value',1.0,'Position',[500 20 100 20],'Callback',@change_factor);
a.uicontrols.reset = uicontrol('style','pushbutton','Callback',@ResetCallback, 'string','RESET','Position',[620 20 100 20]);
a.factor.son_sigma = 0.25;
a.factor.soff_sigma = 0.25;


% UICONTROL for evaluation
a.uicontrol.evaluate = uicontrol('style','pushbutton','Callback',@EvaluateCallback, 'string','EVALUATE','Position',[750 20 100 20]);
a.axeshandles.SOFF_resp = axes('position',[400 420 500 125]);
set(a.axeshandles.SOFF_resp,'XTick',[],'YTick',[],'Box','on'); 
a.axeshandles.SON_resp = axes('position',[400 590 500 125]);
set(a.axeshandles.SON_resp,'XTick',[],'YTick',[],'Box','on');
set(gcf,'UserData',a);
drawnow;

show_stim1();
show_stim2();
show_soffrf();
show_sonrf();

end

function show_stim1(~,~)
a = get(gcf,'Userdata');
a.template.costerm = cos(2*pi*a.template.Y*get(a.uicontrols.sf1,'value')+get(a.uicontrols.phase1,'value'));
img = 0.5*(a.template.costerm.*a.template.expterm);
stim_lms = [get(a.uicontrols.L1,'value'), get(a.uicontrols.M1,'value'), get(a.uicontrols.S1,'value')];
set(a.uicontrols.textL1,'String', strcat('L:',num2str(stim_lms(1))));
set(a.uicontrols.textM1,'String', strcat('M:',num2str(stim_lms(2))));
set(a.uicontrols.textS1,'String', strcat('S:',num2str(stim_lms(3))));
set(a.uicontrols.textsf1,'String', strcat('sf:',num2str(get(a.uicontrols.sf1,'value'))));
set(a.uicontrols.textphase1,'String', strcat('phase:',num2str(get(a.uicontrols.phase1,'value'))));
stim_rgb = stim_lms * a.param.invM; % converting LMS light to RGB light
A = cat(3,stim_rgb(1)*img,stim_rgb(2)*img,stim_rgb(3)*img); a.stim1 = A;
A = (0.5*A./(max(abs(A(:)))+0.01)) + 0.5;
image(A,'Parent',a.axeshandles.stim1);
set(a.axeshandles.stim1,'XTick',[],'YTick',[],'Box','on');
set(gcf,'UserData',a);
end

function show_stim2(~,~)
a = get(gcf,'Userdata');
a.template.costerm = cos(2*pi*a.template.Y*get(a.uicontrols.sf2,'value')+get(a.uicontrols.phase2,'value'));
img = 0.5*(a.template.costerm.*a.template.expterm);
stim_lms = [get(a.uicontrols.L2,'value'), get(a.uicontrols.M2,'value'), get(a.uicontrols.S2,'value')];
set(a.uicontrols.textL2,'String', strcat('L:',num2str(stim_lms(1))));
set(a.uicontrols.textM2,'String', strcat('M:',num2str(stim_lms(2))));
set(a.uicontrols.textS2,'String', strcat('S:',num2str(stim_lms(3))));
set(a.uicontrols.textsf2,'String', strcat('sf:',num2str(get(a.uicontrols.sf2,'value'))));
set(a.uicontrols.textphase2,'String', strcat('phase:',num2str(get(a.uicontrols.phase2,'value'))));
stim_rgb = stim_lms * a.param.invM; % converting LMS light to RGB light
A = cat(3,stim_rgb(1)*img,stim_rgb(2)*img,stim_rgb(3)*img); a.stim2 = A;
A = (0.5*A./(max(abs(A(:)))+0.01)) + 0.5;
image(A,'Parent',a.axeshandles.stim2);
set(a.axeshandles.stim2,'XTick',[],'YTick',[],'Box','on');
set(gcf,'UserData',a);
end

function show_soffrf(~,~)
a = get(gcf,'UserData');
sigma = get(a.uicontrols.sigma_soff,'value');
set(a.uicontrols.textsigma_soff,'String', strcat('sigma:',num2str(sigma)));
a.factor.soff_sigma = sigma;
stim_lms = [0.23, -0.17, -0.32]; % From Tailby et al; 2008
stim_rgb = stim_lms * a.param.invM;
img = exp(-(a.template.X.^2 + a.template.Y.^2)/2/sigma^2);
A = cat(3,stim_rgb(1)*img,stim_rgb(2)*img,stim_rgb(3)*img);
A = (0.5*A./(max(abs(A(:)))+0.01)) + 0.5;
a.SOFF_rgb = stim_rgb;
image(A,'Parent',a.axeshandles.SOFF);
set(a.axeshandles.SOFF,'XTick',[],'YTick',[],'Box','on');
set(gcf,'UserData',a);
end

function show_sonrf(~,~)
a = get(gcf,'UserData');
sigma = get(a.uicontrols.sigma_son,'value');
set(a.uicontrols.textsigma_son,'String', strcat('sigma:',num2str(sigma)));
a.factor.son_sigma = sigma;
stim_lms = [-0.29, -0.05, 0.56]; % From Tailby et al; 2008
stim_rgb = stim_lms * a.param.invM;
img = exp(-(a.template.X.^2 + a.template.Y.^2)/2/sigma^2);
A = cat(3,stim_rgb(1)*img,stim_rgb(2)*img,stim_rgb(3)*img);
A = (0.5*A./(max(abs(A(:)))+0.01)) + 0.5;
a.SON_rgb = stim_rgb;
image(A,'Parent',a.axeshandles.SON);
set(a.axeshandles.SON,'XTick',[],'YTick',[],'Box','on');
set(gcf,'UserData',a);
end

function change_factor(~,~)
a = get(gcf,'UserData');
factor = get(a.uicontrols.factor,'value');
set(a.uicontrols.text_factor,'String', strcat('factor:',num2str(factor)));

% S-OFF
sigma1 = a.factor.soff_sigma;
set(a.uicontrols.sigma_soff,'value',sigma1*factor);
show_soffrf();

% S-ON
sigma2 = a.factor.son_sigma;
set(a.uicontrols.sigma_son,'value',sigma2*factor);
show_sonrf();

set(gcf,'UserData',a);
end

function ResetCallback(~,~)
a = get(gcf,'UserData');

% factor
set(a.uicontrols.text_factor,'String', strcat('factor:',num2str(1)));
set(a.uicontrols.factor, 'value', 1);

% S-OFF
a.factor.soff_sigma = 0.25;
set(a.uicontrols.sigma_soff,'value',a.factor.soff_sigma);
show_soffrf();

% S-ON
a.factor.son_sigma = 0.25;
set(a.uicontrols.sigma_son,'value',a.factor.son_sigma);
show_sonrf();

set(gcf,'UserData',a);
end

function EvaluateCallback(~,~)
a = get(gcf,'UserData');
factors = logspace(0.1,1,21);
stim_effect_SON = []; stim_effect_SOFF = [];

for ii = 1:numel(factors)

    % SON template for different factors
    imgSONrf_template = exp(-(a.template.X.^2 + a.template.Y.^2)/2/(a.factor.son_sigma*factors(ii))^2);
    imgSONrf = cat(3,a.SON_rgb(1)*imgSONrf_template,a.SON_rgb(2)*imgSONrf_template,a.SON_rgb(3)*imgSONrf_template);
    
    % SOFF template varying the factors
    imgSOFFrf_template = exp(-(a.template.X.^2 + a.template.Y.^2)/2/(a.factor.soff_sigma*factors(ii))^2);
    imgSOFFrf = cat(3,a.SOFF_rgb(1)*imgSOFFrf_template,a.SOFF_rgb(2)*imgSOFFrf_template,a.SOFF_rgb(3)*imgSOFFrf_template);
    
    % calculating the dot product
    stim_effect_SON = [stim_effect_SON; dot(a.stim1(:),imgSONrf(:)) dot(a.stim2(:),imgSONrf(:))];
    stim_effect_SOFF = [stim_effect_SOFF; dot(a.stim1(:),imgSOFFrf(:)) dot(a.stim2(:),imgSOFFrf(:))];
end

axes(a.axeshandles.SON_resp); 
plot(stim_effect_SON(:,1), '--ro', 'Linewidth',2); hold on;
plot(stim_effect_SON(:,2), '--go', 'Linewidth',2); 
title('S ON response'); legend('stim1','stim2');
set(a.axeshandles.SON_resp,'xtickLabel', log10(factors));hold off;

axes(a.axeshandles.SOFF_resp);  
plot(stim_effect_SOFF(:,1), '--ro', 'Linewidth',2); hold on;
plot(stim_effect_SOFF(:,2), '--go', 'Linewidth',2); 
title('S OFF response');legend('stim1','stim2'); xlabel('RF size -->');
set(a.axeshandles.SOFF_resp,'xtickLabel', log10(factors));hold off;

set(gcf,'UserData',a);
end