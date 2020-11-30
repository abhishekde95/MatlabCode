function GLMSPopGUI_Mean2Var()

figure(121); clf;
set(gcf,'pos',[200 200 800 500],'NumberTitle','off','Name','Population Mean to Variance')

% Partition figure into panels
popm2vfig.controlpanel = uipanel('parent',gcf,'units','normalized',...
    'pos',[.01 .01 .34 .34]);
popm2vfig.cellspanel = uipanel(gcf,'Units','normalized',...
    'position',[.01 .36 .34 .63]);
papm2vfig.disppanel = uipanel(gcf,'units','normalized',...
    'pos',[.36 .01 .63 .98]);

% Set up control panel
conpanel.analyze = uicontrol(popm2vfig.controlpanel,'units','normalized',...
    'pos',[.1 .05 .4 .2],'string','Analyze','callback',@startanalysis);

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])

% Set up display panel
datafiles = GLMSPopData(:,1);
cellspanel.celltable = uitable(popm2vfig.cellspanel,'units','normalized',...
    'pos',[.01 .01 .98 .98],'data',datafiles)


UnpackPopParams()

% Save figure data
set(popm2vfig.controlpanel,'userdata',conpanel);
set(popm2vfig.disppanel,'userdata',disppanel);
set(popm2vfig.cellpanel,'userdata',cellpanel);
set(gcf,'userdata',m2vfig);


end


function UnpackPopParams()


keyboard

end