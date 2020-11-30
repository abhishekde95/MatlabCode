function [diffs, p] = prettypaircomp(X,varnames,type)
% function [diffs, p] = prettypaircomp(X,<varnames>,<type>)
%
% Makes an interactive image of a pairwise mean difference.
% Click on a cell to see a histogram of differences. (and a scatter plot?)
% Modified from 'prettycor.m'.
% GDLH 12/1/16

k = size(X,2);
if (nargin < 3)
    type = 't-test';
end
if (~any(strcmp(type,{'t-test','wilcoxon'})))
    error('type (3rd arg) must be "t-test" (the default), or "wilcoxon"');
end
if (nargin == 1 )
    varnames = cell(k,1);
end
if (length(varnames) ~= k)
    error('mismatch in # of variables and  labels');
end
diffs = nan*ones(size(X,2),size(X,2));
p = nan*ones(size(X,2),size(X,2));
for i = 1:size(X,2)
    for j = 1:size(X,2)
        if strcmp(type,'t-test')
            diffs(i,j) = nanmean(X(:,i)-X(:,j));
            [h,p(i,j)] = ttest(X(:,i)-X(:,j));
        else % type is wilcoxon
            diffs(i,j) = nanmedian(X(:,i)-X(:,j));
            p(i,j) = signrank((X(:,i)-X(:,j)));
        end
    end
    if strcmp(type,'t-test') % statname used by Makeplot and dispscatter, below
        statname = 'mean';
    else
        statname = 'median';
    end
end

Makeplot(diffs,p);

    function Makeplot(r,p)
        figure; axes; hold on;
        normalized_diffs = (r-min(r(:)))/(max(r(:))-min(r(:)));
        him = image(normalized_diffs*64);  colormap(gray(64));
        axis ij;
        for i = 1:k
            for j = 1:k
                if (i == j)
                    continue
                end
                h = text(i,j-.2,[statname,': ',num2str(r(i,j),2)],'HorizontalAlignment','center','ButtonDownFcn',@dispscatter);
                if (normalized_diffs(i,j)>0.5)
                    set(h,'color',[1 1 1]);
                end
                if (p(i,j) < 0.05)
                    set(h,'color',[1 0 0]);
                end
                
                h = text(i,j+.2,['p = ',num2str(p(i,j),2)],'HorizontalAlignment','center','ButtonDownFcn',@dispscatter);
                if (normalized_diffs(i,j)>0.5)
                    set(h,'color',[1 1 1]);
                end
                if (p(i,j) < 0.05)
                    set(h,'color',[1 0 0]);
                end
            end
        end
        set(gca,'XTick',1:k,'XTickLabel',varnames);
        set(gca,'YTick',1:k,'YTickLabel',varnames);
        axis tight
        set(him,'ButtonDownFcn',@dispscatter);
    end

    % Support function. 
    function dispscatter(h,ev)
        whichpt = get(gca,'CurrentPoint');
        whichpt = round(whichpt(1,[1 2]));
        whichpt = min([whichpt; k k]);
        figure('position',[ 78   517   239   440]);
        subplot(2,1,1); hold on;
        plot(X(:,whichpt(1)),X(:,whichpt(2)),'k.');
        xlabel(varnames(whichpt(1)));
        ylabel(varnames(whichpt(2)));
        minval = min([X(:,whichpt(1)); X(:,whichpt(2))]);
        maxval = max([X(:,whichpt(1)); X(:,whichpt(2))]);
        set(gca,'Xlim',[minval maxval],'Ylim',[minval maxval])
        plot([minval maxval],[minval maxval],'k-');
        axis square;
        title(['p = ',num2str(p(whichpt(1),whichpt(2)))]);
        subplot(2,1,2); hold on;
        hist(X(:,whichpt(1))-X(:,whichpt(2)),max([10,floor(size(X,1)/10)]));
        %xlabel([varnames(whichpt(1)),' - ',varnames(whichpt(2))],'fontsize',10);
        xlabel('X-Y');
        title([statname,': ',num2str(diffs(whichpt(1),whichpt(2)))]);
    end
end