function initGUI()
    global gl

    % GUI variables
    fieldnames = {'pLpM','mLmM','pLmM','mLpM'};
    labels = {'+L+M' '-L-M' '+L-M' '-L+M'};
    initBGcol = [.6 .6 .6];
    panelsize = 75;
    panelXpos = linspace(45,825,10);
    panelYpos = linspace(265,15,4);
    
    % Set up GUI Pannel
    close all;
    figure(1);
    set(figure(1),'DefaultAxesUnits','normalized','position',[100 100 1200 800],'UserData',[])
    a = get(figure(1),'UserData');
    
    % Set up subpannels
    a.axes.conpanel = uipanel('Title','Control Panel','Position',[.015 .515 .185 .480],'BorderType','beveledin','BackgroundColor',[.9 .9 .9]);
    a.axes.dnpanel = uipanel('Title','Dense Noise','Position',[.21 .515 .775 .480],'BorderType','beveledin','BackgroundColor',[.9 .9 .9]);
    a.axes.glmppanel = uipanel('Title','GLMP','Position',[.015 .015 .970 .490],'BorderType','beveledin','BackgroundColor',initBGcol);
    
    
    % Set Up Control Panel
    % Subunit Selector
    conpanel = get(a.axes.conpanel,'UserData');
    conpanel.axes.selector = uibuttongroup('Parent',a.axes.conpanel,'Units','normalized',...
        'visible','off','Position',[.05 .75 .5 .2],'Title','Stimulus Selector',...
        'BorderType','etchedout','Visible','on','BackgroundColor',[.9 .9 .9],...
        'Title','Subunit Selection','TitlePosition','centertop');
    conpanel.uicontrols.sub1 = uicontrol('Parent',conpanel.axes.selector,'Units','normalized',...
        'Style','radiobutton','String','Subunit 1','FontSize',13,...
        'pos',[.05 .6 .9 .3],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.sub2 = uicontrol('parent',conpanel.axes.selector,'Units','normalized',...
        'Style','radiobutton','String','Subunit 2','FontSize',13,...
        'pos',[.05 .2 .9 .3],'BackgroundColor',[.9 .9 .9]);
    %set(conpanel.axes.selector,'SelectedObject',conpanel.uicontrols.sub1,'SelectionChangeFcn',@switchstim);
    
    % Stimulus adjustments (lum/col on/off)
    conpanel.axes.stimpar = uipanel('Parent',a.axes.conpanel,'Position',[.6 .75 .35 .2],'Title','Stim Params',...
        'BorderType','etchedout','BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.lumon = uicontrol('parent',conpanel.axes.stimpar,'units','normalized',...
        'position',[.05 .5 .9 .4],'style','pushbutton','string','Lum On','backgroundcolor',[.5 1 .5],...
        'callback',@lumonoff);
    conpanel.uicontrols.colon = uicontrol('parent',conpanel.axes.stimpar,'units','normalized',...
        'position',[.05 .05 .9 .4],'style','pushbutton','string','Col On',...
        'backgroundcolor',[.5 1 .5],'callback',@colonoff);
    conpanel.uicontrols.reset = uicontrol('parent',a.axes.conpanel,'units','normalized',...
        'pos',[.05 .15 .4 .1],'style','pushbutton','string','Reset STA',...
        'backgroundcolor',[.8 .5 .5],'callback',@resetSTA);
    
    % Epoch Switch
    conpanel.uicontrols.startnextepoch = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','pushbutton','Position',[.05 .01 .9 .1],'string','Start GLMP',...
        'backgroundcolor',[.5 1 .5],'Callback',@startGLMP);
    conpanel.uicontrols.addGabor = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','pushbutton','Position',[.5 .15 .45 .1],'string','Add Gabors',...
        'backgroundcolor',[1 .5 .5],'Callback',@addGabors);

    % Counters
    b.spikecounter = 0;
    b.nstimleft = [];
    b.nstimpresented = 0;
    conpanel.uicontrols.spikecounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','0','FontSize',12,...
        'Position',[.05 .54 .425 .05],'BackgroundColor',[.8 .8 .8],...
        'UserData',b.spikecounter);
    conpanel.labels.spikecounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Number of Spikes Collected','FontSize',10,...
        'HorizontalAlignment','center',...
        'Position',[.05 .6 .425 .07],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.stimqueue = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','~','FontSize',12,...
        'Position',[.525 .54 .425 .05],'BackgroundColor',[.8 .8 .8],...
        'Visible','on','UserData',b.nstimleft);
    conpanel.labels.stimqueue = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Number of Stimuli in Queue','FontSize',10,...
        'HorizontalAlignment','center',...
        'Visible','on','Position',[.525 .6 .425 .07],'BackgroundColor',[.9 .9 .9]);
    conpanel.uicontrols.stimcounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','edit','string','0','FontSize',12,...
        'Position',[.05 .37 .425 .05],'BackgroundColor',[.8 .8 .8],...
        'Visible','on','UserData',b.nstimpresented);
    conpanel.labels.stimcounter = uicontrol('Parent',a.axes.conpanel,'Units','normalized',...
        'style','text','string','Number of Frames Presented','FontSize',10,...
        'HorizontalAlignment','center',...
        'Visible','on','Position',[.05 .43 .425 .07],'BackgroundColor',[.9 .9 .9]);
    
    %Save User Data
    set(a.axes.conpanel,'UserData',conpanel)
    
    % Set Up Dense Noise Panel
    dnpanel = get(a.axes.dnpanel,'UserData');
        
    % Set up Dense Noise figures
    for ypos = 1:numel(panelYpos)
        field = fieldnames{ypos};
        for xpos = 1:numel(panelXpos)
            dnpanel.axes.(field)(xpos) = axes('Parent',a.axes.dnpanel,'Units','pixels',...
                'Xtick',[],'Ytick',[],'box','on',...
                'Position',[panelXpos(xpos) panelYpos(ypos) panelsize panelsize]);
            if xpos == 1
                figure(1);
                dnpanel.labels.(field) = text(-.2,.5,(labels{ypos}),'Rotation',90,'Units','normalized',...
                    'HorizontalAlignment','center');
            end
        end
    end
    
    % Save User Data
    set(a.axes.dnpanel,'UserData',dnpanel)
    
    % GLMP Panel
    figure(1);
    glmppanel = get(a.axes.glmppanel,'UserData');
    glmppanel.axes.sub1 = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'OuterPosition',[.05 .05 .3 .9],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 25],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'box','on');
    set(glmppanel.axes.sub1,'Color',[.3 .3 .3])
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate')
    title('Subunit 1','FontSize',12,'FontWeight','bold')
    
    glmppanel.axes.sub2 = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'OuterPosition',[.35 .05 .3 .9],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 25],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'box','on');    
    set(glmppanel.axes.sub2,'Color',[.3 .3 .3])
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate')
    title('Subunit 2','FontSize',12,'FontWeight','bold')
    
    glmppanel.axes.gabors = axes('Parent',a.axes.glmppanel,'Units','normalized',...
        'OuterPosition',[.65 .05 .3 .9],...
        'XLim',[-gl.maxcc gl.maxcc],'YLim',[-gl.maxcc gl.maxcc],'ZLim',[0 25],...
        'XGrid','on','YGrid','on','ZGrid','on',...
        'box','on');
    set(glmppanel.axes.gabors,'Color',[.3 .3 .3])
    xlabel('Lcc'); ylabel('Mcc'); zlabel('Firing Rate')
    title('Gabors','FontSize',12,'FontWeight','bold')

    % Save User Data
    set(a.axes.glmppanel,'UserData',glmppanel)
    set(figure(1),'UserData',a)
    
    function lumonoff(a,~)
        %newcc = get(conpanel.uicontrols.lumslider,'value');
        %set(conpanel.axes.lumdisplay,'string',newcc)
        %set(a.axes.conpanel,'UserData',conpanel)
        if strcmp(get(a,'String'),'Lum On')
            gl.templum = .001;
            set(a,'BackgroundColor',[1 .5 .5]);
            set(a,'String','Lum Off');
        elseif strcmp(get(a,'String'),'Lum Off')
            gl.templum = gl.lumcc;
            set(a,'BackgroundColor',[.5 1 .5]);
            set(a,'String','Lum On');
        else
            disp('Problem with lum on/off')
            keyboard
        end
    end

    function colonoff(a,~)
        if strcmp(get(a,'String'),'Col On')
            gl.tempcol = .001;
            set(a,'BackgroundColor',[1 .5 .5]);
            set(a,'String','Col Off');
        elseif strcmp(get(a,'String'),'Col Off')
            gl.tempcol = gl.colcc;
            set(a,'BackgroundColor',[.5 1 .5]);
            set(a,'String','Col On');
        else
            disp('Problem with col on/off')
            keyboard
        end
    end

    function resetSTA(~,~)
        gl.reset = 1;
        set(conpanel.uicontrols.reset,'backgroundcolor',[.5 .7 .5],'string','Resetting...')
    end

    function startGLMP(~,~)
        a = get(figure(1),'UserData');
        conpanel = get(a.axes.conpanel,'UserData');
        if any(any(gl.pixelmask.subunit1)) || any(any(gl.pixelmask.subunit2))
            
            if (get(conpanel.uicontrols.startnextepoch,'Value')) == 1;
                set(conpanel.uicontrols.startnextepoch,'BackgroundColor',[0 1 0])
                gl.startnextepoch = 1;
                set(conpanel.uicontrols.startnextepoch,'Enable','off')
            else
                set(conpanel.uicontrols.startnextepoch,'BackgroundColor',[0.9 0.9 0.8])
                gl.startnextepcoh = 0;
            end
        else
            fprintf('\n***** Must Select Subunits Before Proceeding to GLMP! *****\n\n')
        end
    end

    function addGabors(a,~)
        
        if strcmp(get(a,'String'),'Add Gabors')
            gl.addGabors = 1;
            set(a,'BackgroundColor',[.5 1 .5]);
            set(a,'String','Gabors Added');
            disp('Adding gabors to the stimulus set...')
        elseif strcmp(get(a,'String'),'Gabors Added')
            gl.addGabors = 0;
            set(a,'BackgroundColor',[1 .5 .5]);
            set(a,'String','Add Gabors');
            disp('Removing gabors from the stimulus set...')
        end
        
    end

end