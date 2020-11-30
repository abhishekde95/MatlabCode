% Quick script to plot the relevent monitor spectra for the DTMacPig paper
for calname = {'Dell4BitsCal.mat' 'ViewPixx.mat'}
    figure;
    cals = load(calname{:});
    cal = cals.cals{end};
    colors = 'rgb';
    for channel = 2:3
        spds = reshape(cal.rawdata.mon(:,channel), cal.S_device(3), []);
        [u,s,~] = svd(spds);
        subplot(1,3,1); hold on
        plot(diag(s(1:4,1:4)), ['o' colors(channel)])
        set(gca, 'xtick', [], 'yscale', 'log', 'ylim', [0 .1])
        subplot(1,3,2); hold on
        plot(linspace(380,780,cal.S_device(3)),u(:,2),['-' colors(channel)]);
        set(gca,'xlim',[400 600],'ylim',[-0.3 0.3],'ytick',[-0.3 0 0.3])
        subplot(1,3,3); hold on
        plot(linspace(380,780,cal.S_device(3)),spds,['-' colors(channel)]);
        set(gca,'xlim',[400 600])
    end
end
