%Called by LMTFpaperFigures - this compares the sensitivity of monkeys and
%humans as in paper figure 9 (normalized contrast sensitivity vs temporal
%frequency, luminance and chromatic plotted separately

function vals = senscomp_mvh()
mode = 5;
RF_ref = [5 0];
angletol = pi/7;%0.5236;
vals = {};
try
    load 'C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Zack\IsoSamp\private\data\LMTF.mat'
catch
    load('/Users/horwitzlab/Desktop/MatlabCode/Zack/IsoSamp/private/data/LMTF.mat')
end
SIDs = {'A','U','E','G'};
colors = {[.37 .24 .6], [.7 .67 .82], [.9 .38 .003], [.99 .72 .39]};
lumtfs = logspace(log10(1),log10(30),100);
rgtfs = logspace(log10(1),log10(20),100);
figure; axes; hold on;
h = []; p = []; q = [];
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
for j = 1:length(SIDs)
    s = eval(SIDs{j});
    globalparams = getfield(s.legacy,['mode',num2str(mode),'params']);
    model = LMTF_global_to_local_model(globalparams, RF_ref(1), RF_ref(2), mode);
    %disp(model(13))
    %plotting the final model fit lines
    %lumf1 = @(omega)(model(1)*abs(((1i*2*pi*10^model(5).*omega+1).^-model(3))-model(2)*((1i*2*pi*10^(model(5)+model(6)).*omega+1).^-(model(3)+model(4)))));
    %rgf1 = @(omega)(model(1+6)*abs(((1i*2*pi*10^model(5+6).*omega+1).^-model(3+6))-model(2+6)*((1i*2*pi*10^(model(5+6)+model(6+6)).*omega+1).^-(model(3+6)+model(4+6)))));
    %h(j) = plot(lumtfs,lumf1(lumtfs),'-','LineWidth',2,'Color',colors{j});
    %plot(rgtfs,rgf1(rgtfs),':','LineWidth',2,'Color',colors{j});
    % Collecting raw data
    tot_rfs_query = sprintf('SELECT DISTINCT rfX, rfY FROM LMTF WHERE subjID = ''%s'' AND quality = 1 AND recDate > ''2016-05-27'' AND recDate < ''2017-02-03''', SIDs{j});
    tot_rfs = fetch(conn, tot_rfs_query);
    X_lums = {}; Y_lums = {}; X_chrs = {}; Y_chrs = {};
    for i = 1:length(tot_rfs)
        RF_loc = [(tot_rfs{i,1})/10, (tot_rfs{i,2})/10];
        L = all(s.eccs == repmat(10*[RF_loc(1) RF_loc(2)],size(s.eccs,1),1),2);
        if sum(L) ~= 0
            raw = s.raw{L};
            Loog = logical(raw(:,4));
            [th,r] = cart2pol(raw(:,1),raw(:,2));
            scale_model = LMTF_global_to_local_model(globalparams, RF_loc(1), RF_loc(2), mode);
            lumf1 = @(omega)(scale_model(1)*abs(((1i*2*pi*10^scale_model(5).*omega+1).^-scale_model(3))-scale_model(2)*((1i*2*pi*10^(scale_model(5)+scale_model(6)).*omega+1).^-(scale_model(3)+scale_model(4)))));
            rgf1 = @(omega)(scale_model(1+6)*abs(((1i*2*pi*10^scale_model(5+6).*omega+1).^-scale_model(3+6))-scale_model(2+6)*((1i*2*pi*10^(scale_model(5+6)+scale_model(6+6)).*omega+1).^-(scale_model(3+6)+scale_model(4+6)))));
            lum_peak = max(lumf1(lumtfs));
            chr_peak = max(rgf1(rgtfs));
            lum_scale_factor = 1/lum_peak; %model(1)/scale_model(1);
            chr_scale_factor = 1/chr_peak;
            % Getting raw data near the luminance axis
            lum_lim = th<pi/4+angletol/2 & th>pi/4-angletol/2;
            lum_lim = lum_lim & ~Loog;
            % Getting raw data near the chromatic axis
            chr_lim = th<3*pi/4+angletol/2 & th>3*pi/4-angletol/2;
            chr_lim = chr_lim & ~Loog;
            if sum(lum_lim) ~= 0
                try
                    X_lums{j} = [X_lums{j}; raw(lum_lim,3)]; %tf
                    Y_lums{j} = [Y_lums{j}; (1./r(lum_lim))*lum_scale_factor]; %contrast sensitivity
                catch
                    X_lums{j} = raw(lum_lim,3);
                    Y_lums{j} = (1./r(lum_lim))*lum_scale_factor;
                end
            end
            if sum(chr_lim) ~= 0
                try
                    X_chrs{j} = [X_chrs{j}; raw(chr_lim,3)];
                    Y_chrs{j} = [Y_chrs{j}; 1./r(chr_lim)*chr_scale_factor];
                catch
                    X_chrs{j} = raw(chr_lim,3);
                    Y_chrs{j} = 1./r(chr_lim)*chr_scale_factor;
                end
            end
        end
    end
    %Binning and plotting raw data
    p(j) = subplot(1,2,1); hold on;
    set(p(j), 'Yscale', 'log', 'Xscale', 'log');
    set(p(j),'Xlim',[1 30]);
    set(p(j),'Ylim',[0 5]);
    lumbinedges = logspace(log10(1), log10(30),10);
    deltas = logspace(log10(.03), log10(1), 10);
    lumbincenters = sqrt(lumbinedges(1:end-1).*lumbinedges(2:end));
    X_lums_binned = zeros(length(lumbinedges)-1,1); Y_lums_binned = zeros(length(lumbinedges)-1,1); Y_lums_std = zeros(length(lumbinedges)-1,1);
    for i = 1:length(lumbinedges)-1
        L = X_lums{j} >= lumbinedges(i) & X_lums{j} < lumbinedges(i+1);
        if L == 0
            continue;
        else
            X_lums_binned(i) = lumbincenters(i);
            Y_lums_binned(i) = mean(Y_lums{j}(L));
            Y_lums_std(i) = std(Y_lums{j}(L))./sqrt(sum(L));
            if i == 1
                disp(SIDs{j});
                disp(lumbinedges(1:2));
                vals{j} = Y_lums{j}(L);
            end
            if (length(Y_lums{j}(L))) == 1
                plot(X_lums_binned(i), Y_lums_binned(i), 'o', 'MarkerSize', 2, 'MarkerEdgeColor',colors{j}, 'MarkerFaceColor', colors{j});
            else
                plot(X_lums_binned(i), Y_lums_binned(i), 'o', 'MarkerSize', 5,'Color',colors{j});
                se_min = 10.^(mean(log10(Y_lums_binned(i)))-Y_lums_std(i));
                se_max = 10.^(mean(log10(Y_lums_binned(i)))+Y_lums_std(i));
                plot([X_lums_binned(i) X_lums_binned(i)],[se_min se_max],'-', 'Color', colors{j});
                delta = deltas(i);
                plot([X_lums_binned(i)-delta X_lums_binned(i)+delta], [se_max se_max],'-', 'Color', colors{j});
                plot([X_lums_binned(i)-delta X_lums_binned(i)+delta], [se_min se_min],'-', 'Color', colors{j});
            end
        end
    end
    plot(lumtfs, lumf1(lumtfs)*lum_scale_factor,'-','Color',colors{j});
    q(j) = subplot(1,2,2); hold on;
    set(q(j), 'Yscale', 'log', 'Xscale', 'log');
    set(q(j),'Xlim',[1 20]);
    set(q(j),'Ylim',[0 3]);
    chrbinedges = logspace(log10(1), log10(20), 8);
    chrbincenters = sqrt(chrbinedges(1:end-1).*chrbinedges(2:end));
    X_chrs_binned = zeros(length(chrbinedges)-1,1); Y_chrs_binned = zeros(length(chrbinedges)-1,1); Y_chrs_std = zeros(length(chrbinedges)-1,1);
    for i = 1:length(chrbinedges)-1
        L = X_chrs{j} >= chrbinedges(i) & X_chrs{j} < chrbinedges(i+1);
        if i == 1
            vals{j} = [vals{j}; Y_chrs{j}(L)];
        end
        if L == 0
            continue;
        else
            X_chrs_binned(i) = chrbincenters(i);
            Y_chrs_binned(i) = mean(Y_chrs{j}(L));
            Y_chrs_std(i) = std(Y_chrs{j}(L))/sqrt(sum(L));
            if (length(Y_chrs{j}(L))) == 1
                plot(X_chrs_binned(i), Y_chrs_binned(i), 'o', 'MarkerSize', 2, 'MarkerEdgeColor',colors{j}, 'MarkerFaceColor', colors{j});
            else
                plot(X_chrs_binned(i), Y_chrs_binned(i), 'o', 'MarkerSize', 5,'Color',colors{j});
                se_max = 10.^(mean(log10(Y_chrs_binned(i)))+Y_chrs_std(i));
                se_min = 10.^(mean(log10(Y_chrs_binned(i)))-Y_chrs_std(i));
                plot([X_chrs_binned(i) X_chrs_binned(i)],[se_min se_max],'-', 'Color', colors{j});
                delta = deltas(i);
                plot([X_chrs_binned(i)-delta X_chrs_binned(i)+delta], [se_max se_max],'-', 'Color', colors{j});
                plot([X_chrs_binned(i)-delta X_chrs_binned(i)+delta], [se_min se_min],'-', 'Color', colors{j});
            end
        end
    end
    plot(rgtfs,rgf1(rgtfs)*chr_scale_factor,'-','Color',colors{j});
end
%legend([p(1), q(1)], 'lum data', 'chr data'); %h removed, 'A', 'U', 'E', 'G', removed
set([p q], 'Yticklabels', {'1'});
set([p q], 'Xticklabel', {'1', '10'});
set([get(p(j), 'xlabel'), get(q(j), 'xlabel')], 'String', 'Temporal Frequency (Hz)');
set([get(p(j), 'ylabel')], 'String', 'Normalized Contrast Sensitivity');
set([p q], 'TickDir', 'out');
pbaspect(p(j), [1 1 1]);
pbaspect(q(j), [1 1 1]);
title(p(j), 'L+M'); title(q(j), 'L-M'); %title(r(j), 'mean lum bin1, 10 + std err');
end