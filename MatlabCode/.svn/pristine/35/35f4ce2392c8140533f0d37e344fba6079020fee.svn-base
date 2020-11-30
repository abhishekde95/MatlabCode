% This script is a model of a V1 cell tuned to various directions in the LM
% plane.  We want to know how various sampling schemes may bias the
% results.

clear all
close all

% Set up figure
figure(1000); clf;
set(gcf,'units','pixels','pos',[50 100 900 900],'NumberTitle','off',...
    'Name','Sampling Scheme Model');
modelfig = get(gcf,'UserData');
modelfig.conpanel = uipanel('pos',[.75 .625 .225 .35],'parent',gcf,'title','Control Panel');
modelfig.surfpanel = uipanel('Pos',[.025 .625 .7 .35],'Parent',gcf);
modelfig.fitspanel = uipanel('pos',[.025 .025 .95 .575],'parent',gcf);
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');


% Set up controls for analysis
conpanel.uicontrols.nsamps = uicontrol('style','edit','parent',modelfig.conpanel,...
    'pos',[50 250 100 20],'string',1,'fontsize',10);
conpanel.labels.nsamps = uicontrol('Parent',modelfig.conpanel,'Units','pixels',...
    'style','text','string','Samples per Direction','FontSize',10,...
    'HorizontalAlignment','center','Position',get(conpanel.uicontrols.nsamps,'Position')+[-30 20 60 0]);
conpanel.uicontrols.npresstim = uicontrol('parent',modelfig.conpanel,'units','pixels',...
    'style','edit','pos',[50 200 100 20],'string',5,'fontsize',10);
conpanel.labels.npresstim = uicontrol('parent',modelfig.conpanel,'units','pixels',...
    'style','text','string','Samples per (L,M) Value','fontsize',10,...
    'HorizontalAlignment','center','position',get(conpanel.uicontrols.npresstim,'position') + [-30 20 60 0]);
conpanel.uicontrols.whichsamp = uicontrol('parent',modelfig.conpanel,'units','pixels',...
    'style','popupmenu','string',1:str2double(get(conpanel.uicontrols.nsamps,'string')),...
    'pos',[50 150 100 20],'callback',@PlotSample);
conpanel.labels.whichsamp = uicontrol('parent',modelfig.conpanel,'units','pixels',...
    'style','text','string','Examine Sample #:','fontsize',10,...
    'pos',get(conpanel.uicontrols.whichsamp,'position') + [-30 20 60 0]);
conpanel.uicontrols.startanalysis = uicontrol('style','pushbutton','string','Start Analysis',...
    'parent',modelfig.conpanel,'pos',[50 20 100 40],'fontsize',10,'callback',@SSM_StartAnalysis);

% Set up surfpanel
surfpanel.axes.surffig = axes('parent',modelfig.surfpanel,'units','pixels','pos',[50 40 230 230]);
surfpanel.axes.projfig = axes('parent',modelfig.surfpanel,'units','pixels','pos',[350 40 230 230]);

% Set up fitspanel
fitspanel.axes.GOF(1) = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .8 .85 .125],'box','on','xlim',[-pi pi]);
fitspanel.axes.GOF(2) = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .55 .85 .125],'box','on','xlim',[-pi pi]);
fitspanel.axes.GOF(3) = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .3 .85 .125],'box','on','xlim',[-pi pi]);
fitspanel.axes.GOF(4) = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .05 .85 .125],'box','on','xlim',[0 180]);


%Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);
    

