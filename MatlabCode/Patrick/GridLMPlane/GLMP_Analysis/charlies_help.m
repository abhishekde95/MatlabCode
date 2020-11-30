%%
fin

clr = colormap('jet');
rates = normrnd(10, 4, 10,10);
maxRate = max(rates(:));
minRate = min(rates(:));
mapping = linspace(minRate, maxRate, size(clr,1));


figure, hold on
for r = 1:size(rates,1)
    for c = 1:size(rates,2)
        pt = rates(r,c);
        errs = abs(pt-mapping);
        [~,ind] = min(errs);
        plot(r,c,'o','markersize', 15, 'color', clr(ind,:), 'markerfacecolor', clr(ind,:))
    end
end
