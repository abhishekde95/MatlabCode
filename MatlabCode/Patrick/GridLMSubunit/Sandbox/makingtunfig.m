figure(15); clf;
set(15,'pos',[50 75 1200 750],'numbertitle','off','name','Population Tuning');

% Set up panels
TunFig.conpanel = uipanel(15,'units','normalized','pos',[.01 .01 .3 .98]);
TunFig.popanel = uipanel(15,'units','normalized','pos',[.32 .01 .22 .98]);
TunFig.dispanel = uipanel(15,'units','normalized','pos',[.55 .01 .22 .98]);
TunFig.statspanel = uipanel(15,'units','normalized','pos',[.78 .01 .21 .98]);

%%% Set up conpanel %%%
conpanel.uicontrols.reanalyzeall = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Reanalyze All Cells','units','normalized',...
    'pos',[.05 .01 .425 .05],'Foregroundcolor',[1 0 0],'callback',@ReanalAll);
conpanel.uicontrols.reanalyze = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Reanalyze Cell','units','normalized',...
    'pos',[.525 .01 .425 .05],'callback',@Reanal);
conpanel.uicontrols.overview = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Overview','units','normalized',...
    'pos',[.05 .07 .425 .05],'callback',@LoadOverview);
conpanel.uicontrols.analindiv = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Surface Analysis','units','normalized',...
    'pos',[.525 .07 .425 .05],'callback',@callanalysis);
conpanel.uicontrols.analGUI = uicontrol('parent',TunFig.conpanel,...
    'style','pushbutton','string','Load Control Panel','units','normalized',...
    'pos',[.05 .13 .425 .05],'callback',@LoadControlPanel);
strs = {'All Data' 'Nut Data' 'Maui Data'};
conpanel.uicontrols.whichMonk = uicontrol('parent',TunFig.conpanel,...
    'style','popup','units','normalized','pos',[.1 .18 .35 .05],...
    'string',strs,'callback',@ChooseMonk);
strs = {'P-Value' 'Norm LL'};
conpanel.uicontrols.whichMetric = uicontrol(TunFig.conpanel,...
    'style','popup','units','normalized','pos',[.55 .18 .35 .05],...
    'string',strs,'value',2,'callback',@ChooseMetric);
conpanel.selectedidx = [];

conpanel.axes.normLLhist = axes('parent',TunFig.conpanel,'units','normalized',...
    'pos',[.15 .45 .75 .1],'box','on');
conpanel.labels.normLL2D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.05 .38 .25 .03],'string','Norm LL 2D','fontsize',10);
conpanel.labels.normLL1D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.35 .38 .25 .03],'string','Norm LL 1D','fontsize',10);
conpanel.labels.normLLdiff = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.7 .38 .25 .03],'string','Diff Norm LL','fontsize',10);
conpanel.uicontrols.normLL2D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.05 .35 .25 .03]);
conpanel.uicontrols.normLL1D = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.35 .35 .25 .03]);
conpanel.uicontrols.normLLdiff = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.7 .35 .25 .03]);
conpanel.normLLthresh = .1;
conpanel.labels.normLLthresh = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','text','pos',[.2 .3 .4 .03],'string','Norm LL Thresh','fontsize',12);
conpanel.uicontrols.normLLthresh = uicontrol(TunFig.conpanel,'units','normalized',...
    'style','edit','pos',[.6 .3 .25 .03],'string',conpanel.normLLthresh,...
    'Callback',@SelectNLLthresh_setval);


%%% Set up popanel %%%
popanel.axes.oneD = axes('parent',TunFig.popanel,'units','normalized',...
    'pos',[.05 .625 .9 .35]); box on;
popanel.axes.twoD = axes('parent',TunFig.popanel,'units','normalized',...
    'pos',[.05 .125 .9 .35]); box on;
strs = {'All' 'Unidirectional Excitation' 'Bidirectional Excitation' 'Suppressive'};
popanel.uicontrols.which1Dsurf = uicontrol(TunFig.popanel,...
    'style','popup','units','normalized','pos',[.05 .55 .9 .02],...
    'string',strs,'value',1,'callback',@ChooseSurfType);
strs = {'All' 'Pan Color' 'Horseshoe' 'Hypertuned'};
popanel.uicontrols.which2Dsurf = uicontrol(TunFig.popanel,...
    'style','popup','units','normalized','pos',[.05 .05 .9 .02],...
    'string',strs,'value',1,'callback',@ChooseSurfType);
popanel.nangs = 20;

%%% Set up dispanel %%%
dispanel.axes.surface1d = axes('parent',TunFig.dispanel,'units','normalized',...
    'pos',[.2 .7 .7 .25],'box','on'); axis square;
dispanel.axes.surface2d = axes('parent',TunFig.dispanel,'units','normalized',...
    'pos',[.2 .2 .7 .25],'box','on'); axis square;

% 1D params
dispanel.labels.params.A = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.05 .62 .2 .03],'string','A','fontsize',12);
dispanel.uicontrols.params.A = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.25 .62 .2 .03],'fontsize',12);
dispanel.labels.params.bl = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.5 .62 .2 .03],'string','Bl','fontsize',12);
dispanel.uicontrols.params.bl = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.75 .62 .2 .03],'fontsize',12);
dispanel.labels.params.sig = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.05 .58 .2 .03],'string','Sigmas','fontsize',12);
dispanel.uicontrols.params.sig1 = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.25 .58 .2 .03],'fontsize',12);
dispanel.uicontrols.params.sig2 = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.5 .58 .2 .03],'fontsize',12);
dispanel.labels.params.exp = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.25 .51 .2 .03],'string','Exp','fontsize',12);
dispanel.uicontrols.params.exp = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.25 .54 .2 .03],'fontsize',12);
dispanel.labels.params.rot = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.5 .51 .2 .03],'string','Rot','fontsize',12);
dispanel.uicontrols.params.rot = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.5 .54 .2 .03],'fontsize',12);
dispanel.labels.params.kappa = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.75 .51 .2 .03],'string','Kappa','fontsize',12);
dispanel.uicontrols.params.kappa = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.75 .54 .2 .03],'fontsize',12);

% 2D params
dispanel.labels.params.A = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.05 .12 .2 .03],'string','A','fontsize',12);
dispanel.uicontrols.params.A = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.25 .12 .2 .03],'fontsize',12);
dispanel.labels.params.bl = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.5 .12 .2 .03],'string','Bl','fontsize',12);
dispanel.uicontrols.params.bl = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.75 .12 .2 .03],'fontsize',12);
dispanel.labels.params.sig = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.05 .08 .2 .03],'string','Sigmas','fontsize',12);
dispanel.uicontrols.params.sig1 = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.25 .08 .2 .03],'fontsize',12);
dispanel.uicontrols.params.sig2 = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.5 .08 .2 .03],'fontsize',12);
dispanel.uicontrols.params.orthsig = uicontrol(TunFig.dispanel,'units','normalized',...
   'style','edit','pos',[.75 .08 .2 .03],'fontsize',12);
dispanel.labels.params.exp = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.25 .01 .2 .02],'string','Exp','fontsize',12);
dispanel.uicontrols.params.exp = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.25 .04 .2 .03],'fontsize',12);
dispanel.labels.params.rot = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.5 .01 .2 .02],'string','Rot','fontsize',12);
dispanel.uicontrols.params.rot = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.5 .04 .2 .03],'fontsize',12);
dispanel.labels.params.kappa = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','text','pos',[.75 .01 .2 .02],'string','Kappa','fontsize',12);
dispanel.uicontrols.params.kappa = uicontrol(TunFig.dispanel,'units','normalized',...
    'style','edit','pos',[.75 .04 .2 .03],'fontsize',12);


%%% Set up stats panel %%%




