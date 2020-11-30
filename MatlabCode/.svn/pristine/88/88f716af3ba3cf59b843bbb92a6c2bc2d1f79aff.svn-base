%%
% LMTF data and model browser
% Here's how to call this function:
% load /Users/greghorwitz/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat
% LMTFBrowser(A)
%
% This function looks for the "models" field of A. If it doesn't exist,
% but there is a "legacy.mode3models" field, it tries to use this.
% This architecture  allows the user to create a "models" field anyway they want (e.g.
% A.legacy.model1models from LMTF_generate_module_data.m) and use this
% browser to look at residuals from those models.

function LMTFBrowser(LMTFstruct)
UNSELECTEDCOLOR = [0 0 0];
MAXTF = 60;
MINTF = 1;
figure('Units','inches','position',[1 1 11 9]);
% Axes 1: The representation of visual space
hax(1) = axes('units','inches','position',[6 4.5 3 3]); hold on;
spothandles = [];
for i = 1:size(LMTFstruct.eccs,1)
    spothandles(i)=plot(LMTFstruct.eccs(i,1),LMTFstruct.eccs(i,2),'ko','MarkerSize',20,'MarkerFaceColor',UNSELECTEDCOLOR);
    set(spothandles(i),'ButtonDownFcn',{@DisplayTCSCallBack i spothandles(i)});
end
% Need to create "models" field if it is not already present.
% Using mode 5 (yoked tilted rampy trough) by default
if ~isfield(LMTFstruct,'models')
    if isfield(LMTFstruct.legacy,'mode5models')
        LMTFstruct.models = LMTFstruct.legacy.mode5models;
    end
end

axis(hax(1),'equal');
% Axes 2:  % Where the TCS curve will appear
hax(2) = axes('units','inches','position',[2 1 3 3],'Yscale','log'); hold on;
ylabel('sensitivity'); xlabel('tf (Hz)');
% Axes 3:  % Number of data points per retinal location (small, non-
% interactive plot).
hax(3) = axes('units','inches','position',[9.25 5.5 1.5 1.5],'XTick',[],'Ytick',[]); hold on; % Number of data points per retinal location
ns = [];
for i = 1:size(LMTFstruct.eccs,1)
    plot(LMTFstruct.eccs(i,1),LMTFstruct.eccs(i,2),'ko','MarkerSize',numel(LMTFstruct.raw{i})/50,'MarkerFaceColor',[.5 .5 .5]);
end
% Axes 4: Residuals
hax(4) = axes('units','inches','position',[6 1 3 3],'View',[0 90],'Zscale','log','PlotBoxAspectRatio',[1 1 1]);
hold on;

%---------------------------------
% Setting up uicontrols
hcontrols.clearbutton = uibutton('Style','pushbutton','String','clear','Units','inches','Position',[6 7.5 .5 .5],'Callback',@clearfn);
UB = [500 1 50  6 log10(.1)       .5        500 1 50 6 log10(.1)      .5          90];
LB = [0   0  1 0 log10(.001) log10(1.00001) 0   0  1 0 log10(.001) log10(1.00001) 0]; 

sliderlabels = {'xi','zeta','n1','delta n','tau1','kappa','xi','zeta','n1','delta n','tau1','kappa','theta'};
SLIDERHEIGHT = .3; SLIDERMARGIN = 0; SLIDERY = 8; SLIDERLENGTH = 3;
for i = 1:length(sliderlabels)
    hcontrols.paramslider(i) = uicontrol('Style','slider','Min',LB(i),'Max',UB(i),'Value',LB(i),'Units','inches','position',[1 SLIDERY SLIDERLENGTH SLIDERHEIGHT],'Callback',{@SliderCallback i});
    hcontrols.editbox(i) = uicontrol('Style','edit','String',num2str(LB(i)),'Units','inches','position',[SLIDERLENGTH+1 SLIDERY 1 SLIDERHEIGHT],'Callback',{@EditboxCallback i});
    hcontrols.textbox(i) = uicontrol('Style','text','String',sliderlabels{i},'HorizontalAlignment','left','Units','inches','position',[.5 SLIDERY .5 SLIDERHEIGHT]);
    if (i>6)
        set(hcontrols.textbox(i),'ForegroundColor','red');
    end
    SLIDERY=SLIDERY-SLIDERHEIGHT-SLIDERMARGIN;
end


% Nested functions
    function DisplayTCSCallBack(~,~,which,spothandle)
        cmap = hsv(8);
        if any(get(spothandle,'MarkerFaceColor') ~= UNSELECTEDCOLOR) % If we're deselecting
            set(spothandle,'MarkerFaceColor',UNSELECTEDCOLOR)
            delete(get(spothandle,'UserData')); % using the userdata field to hold on to the handles of the sensitivity curves
        else % we're selecting
            spotcolors = cell2mat(get(spothandles,'MarkerFaceColor'));
            nselected = sum(~all(spotcolors == repmat(UNSELECTEDCOLOR,size(spotcolors,1),1),2));
            model = LMTFstruct.models(:,which);
            % Function handles for the two filters
            [f1,f2] = getTCSfunctionhandles(model);
            SELECTEDCOLOR = cmap(min(nselected+1,8),:);
            h = DisplayTCS(f1,f2,SELECTEDCOLOR);
            
            set(spothandle,'MarkerFaceColor',SELECTEDCOLOR);
            set(spothandle,'UserData',h);
            
            % Adjusting the uicontrols
            for i=1:length(sliderlabels)
                try
                    set(hcontrols.paramslider(i),'Value',model(i))
                catch
                    set(hcontrols.paramslider(i),'Value',LB(i))
                    warning('Parameter out of bounds');
                end
                if (strcmp(get(hcontrols.textbox(i),'String'),'theta'))
                    set(hcontrols.editbox(i),'String',[num2str(model(i)*180/pi),' (deg)']); % theta in degrees
                else
                    set(hcontrols.editbox(i),'String',num2str(model(i)));
                end
            end
            cla(hax(4));
            PlotResiduals(which);
            % Adding a figure legend to hax(4)
            text('Units','inches','Position',[2.75 2],'Color','blue','String','Overestimate sens.');
            text('Units','inches','Position',[2.75 1.75],'Color','red','String','Underestimate sens.');
        end
    end

    function [f1, f2] = getTCSfunctionhandles(model)
        if length(model) < 12
            error('Need at least 12 model parameters');
        end
 %       f1 = @(omega)model(1)*abs(((1i*2*pi*model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*model(5)*model(6).*omega+1).^-model(4)));
 %       f2 = @(omega)model(7)*abs(((1i*2*pi*model(11).*omega+1).^-model(9))-model(8)*((1i*2*pi*model(11)*model(12).*omega+1).^-model(10)));
         f1 = @(omega)model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4))));
         f2 = @(omega)model(7)*abs(((1i*2*pi*10^model(11).*omega+1).^-model(9))-model(8)*((1i*2*pi*10^(model(11)+model(12)).*omega+1).^-(model(9)+model(10))));
    end

    function h = DisplayTCS(f1,f2,SELECTEDCOLOR)
        tf = logspace(log10(MINTF),log10(MAXTF),50);
        sensLUM = f1(tf);
        sensRG = f2(tf);
        axes(hax(2));
        h(1) = plot(tf,sensLUM,'k-','LineWidth',2,'Color',SELECTEDCOLOR);
        h(2) = plot(tf,sensRG,'r--','LineWidth',2,'Color',SELECTEDCOLOR);
    end

    function clearfn(a,b)
        cla(hax(2));
        cla(hax(4));
        set(spothandles(:),'MarkerFaceColor',UNSELECTEDCOLOR,'UserData',[]);
    end

    function PlotResiduals(which)
        residsPlot(LMTFstruct.raw{which}(:,[1:4]),LMTFstruct.models(:,which),[],hax(4));
        set(gca,'Xtick',[-1 -.5 -.2 0 .2 .5 1],'Ytick',[-1 -.5 -.2 0 .2 .5 1])
    end

    % ------------
    % UI callbacks
    % ------------
    function SliderCallback(a,b,which)
        num = get(hcontrols.paramslider(which),'Value');
        set(hcontrols.editbox(which),'String',num2str(num));
        UpdateTCSPlot;
    end
    function EditboxCallback(a,b,which)
        num = str2num(get(hcontrols.editbox(which),'String'));
        set(hcontrols.paramslider(which),'Value',num);
        UpdateTCSPlot;
    end
    function UpdateTCSPlot
        c = get(spothandles(:),'MarkerFaceColor');
        Lselectedspots = zeros(length(c),1);
        for i = 1:length(c)
            Lselectedspots(i) = logical(~all(c{i}==UNSELECTEDCOLOR));
        end
        if (sum(Lselectedspots) == 1)
            whichspot = find(Lselectedspots);
            modelparams = [];
            for j = 1:length(hcontrols.paramslider)
                modelparams(j) = get(hcontrols.paramslider(j),'Value');
            end
            [f1,f2] = getTCSfunctionhandles(modelparams);
            delete(get(spothandles(whichspot),'UserData'));
            h = DisplayTCS(f1,f2, c{whichspot});
            set(spothandles(whichspot),'UserData',h);
        end
    end
end