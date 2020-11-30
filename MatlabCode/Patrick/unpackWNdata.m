% for testing the unpacking proceedures of WN data

global DN

stim = {'pLpM'  'mLmM' 'pLmM' 'mLpM'};
nframesback = 10;

panelXpos = linspace(35,750,nframesback);
panelYpos = linspace(250,nframesback,numel(stim));
panelsize = 75;
nstix = DN.NStixGrid(1);

figure(15); clf;
set(gcf,'units','pixels')
% Plot DN data
for s = 1:numel(stim)
    im = reshape(DN.stats.(stim{s}).STA(1:nstix^2,:),[nstix nstix 1 10]);
    im = abs(im(:,:,1,:));
    im = (im - min(im(:)))./(max(im(:))-min(im(:)));
    for t = 1:nframesback
        dnpanel.axes(s,t) = axes('Parent',gcf,'Units','pixels',...
            'Xtick',[],'Ytick',[],'box','on',...
            'Position',[panelXpos(t) panelYpos(s) panelsize panelsize]);
        tempim = im(:,:,:,t);
        tempim = repmat(tempim,1,1,3,1);
        h = image(tempim,'parent',dnpanel.axes(s,t));
    end
end

%%

% Plot White Noise STA over frames
figure(16); clf;
set(gcf,'units','pixels')
for n = 1:numel(stim)
    im = abs(DN.stats.(stim{n}).STA(1:DN.NStixGrid(1)^2,:));
    im = reshape(im,DN.NStixGrid(1),DN.NStixGrid(1),nframesback);
    im = (im - min(im(:)))./(max(im(:))-min(im(:)));
    for t = 1:nframesback
        dnpanel.axes(n,t) = axes('Parent',gcf,'Units','pixels',...
            'Xtick',[],'Ytick',[],'box','on',...
            'Position',[panelXpos(t) panelYpos(n) panelsize panelsize]);
        tempim = im(:,:,t);
        tempim = repmat(tempim,1,1,3);
        h = image(tempim,'parent',dnpanel.axes(n,t));
    end
end
