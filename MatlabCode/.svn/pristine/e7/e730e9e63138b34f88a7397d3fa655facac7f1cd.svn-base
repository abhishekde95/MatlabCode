% Analyzing three pulse Monte data
stro = nex2stro
ntrials = size(stro.trial,1);
intervals = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'pulse_offset'));
correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'correct'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'fpoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'targon_t'));

uniqueintervals = unique(intervals);
for i = 1:length(uniqueintervals)
    L = intervals == uniqueintervals(i);
    [uniqueintervals(i) sum(correct(L))./sum(L)]
end

astart_t = cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells, 'anlgStartTime')));
samplerate = stro.sum.analog.storeRates{1};
delta_t = 1./samplerate;
e1h = stro.ras(:,strcmp(stro.sum.rasterCells, 'AD11'));
e1v = stro.ras(:,strcmp(stro.sum.rasterCells, 'AD12'));

% Analog channels
for WHICHCOL = 1:3
    figure; axes; hold on;
    for i = 1:ntrials
        ntimesteps = length(stro.ras{i,WHICHCOL});
        t = astart_t(i):delta_t:(astart_t(i)+delta_t*(ntimesteps-1));
        t = t-targon_t(i);
        plot(t, stro.ras{i,WHICHCOL});
    end
end