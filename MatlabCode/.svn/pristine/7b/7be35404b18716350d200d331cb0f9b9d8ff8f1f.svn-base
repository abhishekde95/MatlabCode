%% Yasmine El-Shamayleh (4/2015)

% Hunting for Expression in Cerebellum
% AAV1-SL7-ChR2-mCherry
% AAV9-SL7-ChR2-mCherry
% AAV5-hSyn-ArchT-eYFP

%% INJECTIONS

nInj = 3;
inj_coords(1,:) = [-4,-4];   % AAV1-SL7-ChR2-mCherry
inj_coords(2,:) = [+3.5,-4]; % AAV9-SL7-ChR2-mCherry
inj_coords(3,:) = [0, -5];   % AAV5-hSyn-ArchT-eYFP

%% TRACKS


trk_coords(1,:) = [0,-5];
trk_coords(2,:) = [0,-5];
trk_coords(3,:) = [0,-5];
trk_coords(4,:) = [-1,-5];
trk_coords(5,:) = [-0.5,-5];

trk_coords(6,:) = [-4,-4];
trk_coords(7,:) = [-4,-4];
trk_coords(8,:) = [-3, -4];

trk_coords(9,:) = [3.75,-4];
trk_coords(10,:) = [3.75,-4];
trk_coords(11,:) = [3.75,-4];


%% PLOT

figure;

plot(inj_coords(1,1), inj_coords(1,2),'r.','MarkerSize',50); hold on; box off; axis square;
plot(inj_coords(2,1), inj_coords(2,2),'g.','MarkerSize',50);
plot(inj_coords(3,1), inj_coords(3,2),'m.','MarkerSize',50);
set(gca,'xlim',[-7 7], 'ylim', [-7 7]);
set(gca,'xtick',-7:7,'ytick',-7:7,'tickDir','out');


plot([0 0], [-7 7],'--','color',[.5,.5,.5]); hold on; box off; axis square;
plot([-7 7], [0 0],'--','color',[.5,.5,.5]); hold on; box off; axis square;

plot(trk_coords(1,1), trk_coords(1,2),'kx','MarkerSize',20);
plot(trk_coords(2,1), trk_coords(2,2),'ks','MarkerSize',20);
plot(trk_coords(3,1), trk_coords(3,2),'kv','MarkerSize',20);
plot(trk_coords(4,1), trk_coords(4,2),'k*','MarkerSize',20);
plot(trk_coords(5,1), trk_coords(5,2),'ko','MarkerSize',20);


plot(trk_coords(6,1), trk_coords(6,2),'kx','MarkerSize',20);
plot(trk_coords(7,1), trk_coords(7,2),'ks','MarkerSize',20);
plot(trk_coords(8,1), trk_coords(8,2),'kv','MarkerSize',20);

plot(trk_coords(9,1), trk_coords(9,2),'kx','MarkerSize',20);
plot(trk_coords(10,1), trk_coords(10,2),'ks','MarkerSize',20);
plot(trk_coords(11,1), trk_coords(11,2),'kv','MarkerSize',20);




legend('AAV1-SL7-ChR2-mCherry','AAV9-SL7-ChR2-mCherry','AAV5-hSyn-ArchT-eYFP')
title('sedna''s cerebellum chamber: drive coordinates');
xlabel('Lateral to Medial')
ylabel('Posterior to Anterior')






