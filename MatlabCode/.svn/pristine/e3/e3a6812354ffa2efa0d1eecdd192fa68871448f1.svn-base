function GLMSDAnalysisGUI()
global GLMSD

% Plot figure
figure(1); clf;
set(gcf,'Position',[1150 550 250 250],'toolbar','none',...
    'NumberTitle','off','Name',['GLMSD: ' GLMSD.datafile])
GLMSDGUI = get(gcf,'UserData');

GLMSDGUI.analyses.RespMaps = uicontrol('Parent',gcf,'Units','pixels',...
    'style','pushbutton','string','Response Maps',...
    'Position',[20 180 100 40],'Callback',@GLMSDGUI_RespMaps);
GLMSDGUI.analyses.Surfaces = uicontrol('Parent',gcf,'Units','pixels',...
    'style','pushbutton','string','Surfaces',...
    'Position',[20 120 100 40],'Callback',@GLMSDGUI_Surfaces);

% Save user data
set(gcf,'UserData',GLMSDGUI)

end


