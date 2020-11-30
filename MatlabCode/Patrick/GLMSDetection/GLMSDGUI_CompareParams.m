function GLMSDGUI_CompareParams()
global GLMSD

% NOTE: Having trouble sending variable parameters into fmincon.  This
% might involve more custom programing than its worth.

figure(21); clf;
set(gcf,'units','pixels','Position',[200 200 1000 550],...
    'NumberTitle','off','Name',['Surface Analysis (' GLMSD.datafile ')'])
SurfFig.surfpanel1 = uipanel('Parent',gcf,'Units','normalized',...
    'position',[.025 .45 .3 .525],'title','Model 1 Surface');
SurfFig.surfpanel2 = uipanel('Parent',gcf,'Units','normalized',...
    'position',[.35 .45 .3 .525],'title','Model 2 Surface');
SurfFig.diffpanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.675 .45 .3 .525],'title','Difference Between Surfaces');
SurfFig.conpanel(1) = uipanel('parent',gcf,'units','normalized',...
    'position',[.025 .025 .3 .4],'title','Model 1 Parameters');
SurfFig.conpanel(2) = uipanel('parent',gcf,'units','normalized',...
    'position',[.35 .025 .3 .4],'title','Model 2 Parameters');
SurfFig.statspanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.675 .025 .3 .4],'title','Statistics');

% Build control panel
for n = 1:2
    
    % build seperate axis panels
    cpanel = get(SurfFig.conpanel(n),'UserData');
    cpanel.ax1 = uipanel('parent',SurfFig.conpanel(n),...
        'Position',[.025 .5 .95 .45],'title','Axis 1');
    cpanel.ax2 = uipanel('parent',SurfFig.conpanel(n),...
        'position',[.025 .025 .95 .45],'title','Axis 2');
    
    % Ax 1
    % Create ax1 sigma buttons
    cpanel.buttons.ax1.sigpanel = uibuttongroup('parent',cpanel.ax1,...
        'position',[0 .025 .3 .9],'title','Sigma','borderwidth',0,...
        'SelectionChangeFcn',@ChangeParams);
    cpanel.buttons.ax1.sig1 = uicontrol('Style','Radio','String','Sym',...
        'pos',[10 30 100 30],'parent',cpanel.buttons.ax1.sigpanel,'HandleVisibility','off');
    cpanel.buttons.ax1.sig2 = uicontrol('Style','Radio','String','Asym',...
        'pos',[10 10 100 30],'parent',cpanel.buttons.ax1.sigpanel,'HandleVisibility','off');
    
    % Create ax1 exp buttons
    cpanel.buttons.ax1.exppanel = uibuttongroup('parent',cpanel.ax1,...
        'position',[.325 .025 .3 .9],'title','Exponent','borderwidth',0,...
        'SelectionChangeFcn',@ChangeParams);
    cpanel.buttons.ax1.exp1 = uicontrol('Style','Radio','String','Sym',...
        'pos',[10 30 100 30],'parent',cpanel.buttons.ax1.exppanel,'HandleVisibility','off');
    cpanel.buttons.ax1.exp2 = uicontrol('Style','Radio','String','Asym',...
        'pos',[10 10 100 30],'parent',cpanel.buttons.ax1.exppanel,'HandleVisibility','off');
    
    % Create ax1 asymptote buttons
    cpanel.buttons.ax1.asmpanel = uibuttongroup('parent',cpanel.ax1,...
        'position',[.65 .025 .3 .9],'title','Upper Asm','borderwidth',0,...
        'SelectionChangeFcn',@ChangeParams);
    cpanel.buttons.ax1.asm1 = uicontrol('Style','Radio','String','Sym',...
        'pos',[10 30 100 30],'parent',cpanel.buttons.ax1.asmpanel,'HandleVisibility','off');
    cpanel.buttons.ax1.asm2 = uicontrol('Style','Radio','String','Asym',...
        'pos',[10 10 100 30],'parent',cpanel.buttons.ax1.asmpanel,'HandleVisibility','off');
    
    % Ax 2
    cpanel.buttons.ax2.ax2on = uicontrol('style','radio','string','On',...
        'units','normalized','position',[.01 .025 .15 .9],'parent',cpanel.ax2,...
        'Callback',@TurnOnAx2,'value',1);
    
    % Create ax2 sigma buttons
    cpanel.buttons.ax2.sigpanel = uibuttongroup('parent',cpanel.ax2,...
        'position',[.16 .025 .25 .9],'title','Sigma','borderwidth',0,...
        'SelectionChangeFcn',@ChangeParams);
    cpanel.buttons.ax2.sig1 = uicontrol('Style','Radio','String','Sym',...
        'pos',[10 30 100 30],'parent',cpanel.buttons.ax2.sigpanel,'HandleVisibility','off');
    cpanel.buttons.ax2.sig2 = uicontrol('Style','Radio','String','Asym',...
        'pos',[10 10 100 30],'parent',cpanel.buttons.ax2.sigpanel,'HandleVisibility','off');
    
    % Create ax2 exp buttons
    cpanel.buttons.ax2.exppanel = uibuttongroup('parent',cpanel.ax2,...
        'position',[.44 .025 .28 .9],'title','Exponent','borderwidth',0,...
        'SelectionChangeFcn',@ChangeParams);
    cpanel.buttons.ax2.exp1 = uicontrol('Style','Radio','String','Sym',...
        'pos',[10 30 100 30],'parent',cpanel.buttons.ax2.exppanel,'HandleVisibility','off');
    cpanel.buttons.ax2.exp2 = uicontrol('Style','Radio','String','Asym',...
        'pos',[10 10 100 30],'parent',cpanel.buttons.ax2.exppanel,'HandleVisibility','off');
    
    % Create ax2 asymptote buttons
    cpanel.buttons.ax2.asmpanel = uibuttongroup('parent',cpanel.ax2,...
        'position',[.72 .025 .25 .9],'title','Upper Asm','borderwidth',0,...
        'SelectionChangeFcn',@ChangeParams);
    cpanel.buttons.ax2.asm1 = uicontrol('Style','Radio','String','Sym',...
        'pos',[10 30 100 30],'parent',cpanel.buttons.ax2.asmpanel,'HandleVisibility','off');
    cpanel.buttons.ax2.asm2 = uicontrol('Style','Radio','String','Asym',...
        'pos',[10 10 100 30],'parent',cpanel.buttons.ax2.asmpanel,'HandleVisibility','off');
    
    % Set up params
    params.sig = 1;
    params.exp = 1;
    params.asm = 1;
    params.sum = params.sig + params.exp + params.asm;
    
    % Save User Variables
    set(cpanel.ax1,'userdata',params);
    set(cpanel.ax2,'userdata',params);
    set(SurfFig.conpanel(n),'userdata',cpanel)
    
end

% Set up stats panel
spanel.df = uicontrol('style','edit','parent',SurfFig.statspanel,...
    'units','normalized','pos',[.725 .75 .25 .2],'string','0');
spanel.dflabel = uicontrol('style','text','parent',SurfFig.statspanel,...
    'units','normalized','pos',[.025 .75 .675 .15],...
    'string','Difference in Degrees of Freedom =','fontsize',12);
spanel.run = uicontrol('style','pushbutton','parent',SurfFig.statspanel,...
    'units','normalized','pos',[.6 .025 .375 .2],...
    'string','Start Analysis','callback',@StartAnalysis);


% Save figure variables
set(SurfFig.statspanel,'userdata',spanel);
set(gcf,'userdata',SurfFig)

end

function TurnOnAx2(a,~)
    
    SurfFig = get(gcf,'Userdata');
    cpanelh = SurfFig.conpanel(get(get(a,'parent'),'parent')==SurfFig.conpanel);
    cpanel = get(cpanelh,'userdata');
    
    if get(a,'value') == 1
        set(cpanel.buttons.ax2.sig1,'enable','on')
        set(cpanel.buttons.ax2.sig2,'enable','on')
        set(cpanel.buttons.ax2.exp1,'enable','on')
        set(cpanel.buttons.ax2.exp2,'enable','on')
        set(cpanel.buttons.ax2.asm1,'enable','on')
        set(cpanel.buttons.ax2.asm2,'enable','on')
    else
        set(cpanel.buttons.ax2.sig1,'enable','off')
        set(cpanel.buttons.ax2.sig2,'enable','off')
        set(cpanel.buttons.ax2.exp1,'enable','off')
        set(cpanel.buttons.ax2.exp2,'enable','off')
        set(cpanel.buttons.ax2.asm1,'enable','off')
        set(cpanel.buttons.ax2.asm2,'enable','off')
    end
    
    set(cpanelh,'userdata',cpanel)
    set(gcf,'userdata',SurfFig)
    
    % Update df tracker
    RecalcParams()

    
end

function ChangeParams(a,b)

    SurfFig = get(gcf,'Userdata');
    cpanelh = SurfFig.conpanel(get(get(a,'parent'),'parent')==SurfFig.conpanel);
    cpanel = get(cpanelh,'userdata');
    if strcmp(get(get(a,'parent'),'title'),'Axis 1')
        ax = 'ax1';
    else
        ax = 'ax2';
    end
    params = get(cpanel.(ax),'userdata');
    
    if strcmp(get(a,'title'),'Upper Asm')
        if strcmp(get(b.NewValue,'string'),'Sym')
            params.asm = 1;
        else
            params.asm = 2;
        end
    elseif strcmp(get(a,'title'),'Exponent')
        if strcmp(get(b.NewValue,'string'),'Sym')
            params.exp = 1;
        else
            params.exp = 2;
        end  
    elseif strcmp(get(a,'title'),'Sigma')
        if strcmp(get(b.NewValue,'string'),'Sym')
            params.sig = 1;
        else
            params.sig = 2;
        end
    end
    params.sum = params.sig + params.exp + params.asm;
    
    % Save user variables
    set(cpanel.(ax),'userdata',params);
    set(cpanelh,'userdata',cpanel);
    set(gcf,'userdata',SurfFig);
    
    % Update df tracker
    RecalcParams()
    
end

function RecalcParams

    SurfFig = get(gcf,'Userdata');
    cpanel1 = get(SurfFig.conpanel(1),'userdata');
    cpanel2 = get(SurfFig.conpanel(2),'userdata');
    spanel = get(SurfFig.statspanel,'userdata');
    
    % Model 1 params
    params1.ax1 = get(cpanel1.ax1,'userdata');
    params1.ax2 = get(cpanel1.ax2,'userdata');
    
    % Model 2 params
    params2.ax1 = get(cpanel2.ax1,'userdata');
    params2.ax2 = get(cpanel2.ax2,'userdata');
    
    % take difference
    nparams1 = params1.ax1.sum;
    if get(cpanel1.buttons.ax2.ax2on,'value')
        nparams1 = nparams1 + params1.ax2.sum;
    end
    nparams2 = params2.ax1.sum;
    if get(cpanel2.buttons.ax2.ax2on,'value')
        nparams2 = nparams2 + params2.ax2.sum;
    end
    
    set(spanel.df,'string',num2str(abs(nparams1-nparams2)))
    
end


function StartAnalysis(~,~)
global GLMSD

% Load Figure Variables
SurfFig = get(gcf,'UserData');

% Set up x,y,z vals
xvals = GLMSD.Lcc;
yvals = GLMSD.Mcc;
zvals = GLMSD.AnsCorrect;

% Set up some variables
angs = linspace(0,pi,5);
nrots = numel(angs);
GOF = nan(nrots,1);
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');

% Define min and max for parameter range
sigmin = .0001;
sigmax = 5;
expmin = .0001;
expmax = 5;
asmmin = .0001;
asmmax = 1;
blmin = 0;
blmax = 1;
angmin = -pi;
angmax = pi;

% Define guesses for parameters
sigGuess = .025;
expGuess = 2;
asmGuess = 1;
blGuess = .5;

for n = 1:2
    cpanel = get(SurfFig.conpanel(n),'UserData');

    params1 = get(cpanel.ax1,'userdata');
    params2 = get(cpanel.ax2,'userdata');
    
    % Sigma
    sigAx1ID = repmat('sig',1,params1.sig);
    sigAx2ID = repmat('sig',1,params2.sig);
    sigAx1PG = repmat(sigGuess,1,params1.sig);
    sigAx2PG = repmat(sigGuess,1,params2.sig);
    sigAx1LB = repmat(sigmin,1,params1.sig);
    sigAx2LB = repmat(sigmin,1,params2.sig);
    sigAx1UB = repmat(sigmax,1,params1.sig);
    sigAx2UB = repmat(sigmax,1,params2.sig);
    
    % Exp
    expAx1ID = repmat('exp',1,params1.exp);
    expAx2ID = repmat('exp',1,params2.exp);
    expAx1PG = repmat(expGuess,1,params1.exp);
    expAx2PG = repmat(expGuess,1,params2.exp);
    expAx1LB = repmat(expmin,1,params1.exp);
    expAx2LB = repmat(expmin,1,params2.exp);
    expAx1UB = repmat(expmax,1,params1.exp);
    expAx2UB = repmat(expmax,1,params2.exp);
    
    % Asm
    asmAx1ID = repmat('asm',1,params1.asm);
    asmAx2ID = repmat('asm',1,params2.asm);
    asmAx1PG = repmat(asmGuess,1,params1.asm);
    asmAx2PG = repmat(asmGuess,1,params2.asm);
    asmAx1LB = repmat(asmmin,1,params1.asm);
    asmAx2LB = repmat(asmmin,1,params2.asm);
    asmAx1UB = repmat(asmmax,1,params1.asm);
    asmAx2UB = repmat(asmmax,1,params2.asm);
    
    params1 = [sigAx1ID expAx1ID asmAx1ID];
    params2 = [sigAx2ID expAx2ID asmAx2ID];
    axis1 = ones(1,numel(params1));
    axis2 = ones(1,numel(params2))*2;
    params1Guess = [sigAx1PG expAx1PG asmAx1PG];
    params2Guess = [sigAx2PG expAx2PG asmAx2PG];

    params0 = [params1 params2 'bl' 'ang'];
    axis = [axis1 axis2 0 0];
    paramsGuess = [params1Guess params2Guess blGuess 0];
    vlb = [sigAx1LB expAx1LB asmAx1LB sigAx2LB expAx2LB asmAx2LB blmin angmin];
    vub = [sigAx1UB expAx1UB asmAx1UB sigAx2UB expAx2UB asmAx2UB blmax angmax];
    
    % Rotate through angles
    for rot = 1:nrots
        
        % Generating an initial guess
        paramsGuess(end) = angs(rot);
        
        % Fit all rotated points
        [f1,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,'surface7','bernoulli',axis,paramsGuess);
        params{rot} = f1;
        GOF(rot) = fval;
        
    end
    
% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Display surface and pts
x = linspace(-max(xvals),max(xvals),50);
[xx yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface5');
surface = reshape(surface,size(xx));
axes(neuropanel.axes.twoD); cla; hold on;
neuropanel.pts.twoD = plot3(xvals,yvals,zvals,'k*');
neuropanel.surf.twoD = surfc(xx,yy,surface);
alpha(.4); drawnow;
    
% Pull out parameters for display
statstext{1} = rnddec(params1(1) + params1(6),2);
statstext{2} = [num2str(rnddec(params1(2),2)) '  '...
    num2str(rnddec(params1(3),2)) '  '...
    num2str(rnddec(params1(4),2)) '  '...
    num2str(rnddec(params1(5),2))];
statstext{3} = [num2str(rnddec(params1(6),2)) '  '...
    num2str(rnddec(params1(7),2)) '  '...
    num2str(rnddec(params1(8),2)) '  '...
    num2str(rnddec(params1(9),2))];
statstext{4} = rnddec(params1(10),2);
statstext{5} = rnddec(params1(11)/pi*180,2);
statstext{6} = [];
predpts = ComputeNakaRushtonJPW(params1,[xvals yvals],'surface5');
chisq = sum((zvals - predpts).^2);
statstext{7} = chisq;
statstext = statstext';

% Display and save stats
neurostats.twoD.params = params1;
neurostats.twoD.negLL = GOF(bestIdx);
set(neurostats.twoD.paramDisp,'string',statstext);
    
% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(neuropanel.axes.twoD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save user variables
set(SurfFig.disp.neuropanel,'UserData',neuropanel)
set(SurfFig.stats.neuro,'UserData',neurostats)
set(SurfFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',SurfFig)
end


end

