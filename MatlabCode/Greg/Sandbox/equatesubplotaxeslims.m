function equatesubplotaxeslims(equateXY)
% Set all the axes of subplots to be same (picking the largest values)
% Set the optional single argument = 1 if you want the X and Y scales to
% be the same. It defaults to 0.

if nargin == 0
    equateXY = 0;
end

children = get(gcf,'Children');
xlims = [];
ylims = [];
for i = 1:length(children)
    if ~isempty(get(children(i),'Children'))
        xlims = [xlims; get(children(i),'Xlim')];
        ylims = [ylims; get(children(i),'Ylim')];
    end
end

for i = 1:length(children)
    set(children(i),'Xlim',[min(xlims(:,1)) max(xlims(:,2))],'Ylim',[min(ylims(:,1)) max(ylims(:,2))]);
    if equateXY
        set(children(i),'DataAspectRatio',[1 1 1]);
    end
end