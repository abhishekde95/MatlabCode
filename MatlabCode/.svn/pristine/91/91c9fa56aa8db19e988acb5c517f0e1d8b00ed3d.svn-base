% Extracting few parameters from DToneloc opto files: sigma, sf,tf,stim duration & laser duration
% Author - Abhishek De, 5/18
close all; clearvars;
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM DToneloc');
training = fetch(conn,'SELECT training FROM DToneloc');
comments = fetch(conn,'SELECT comments FROM DToneloc');
RWmultiplier = fetch(conn,'SELECT RW_multiplier FROM DToneloc');
laserdial = fetch(conn,'SELECT laser_dial FROM DToneloc');
close(conn);

% Interested in opto files, i.e. training == no
idx  = find(strcmp(training,'no'));
filenameoptoexp = filename(idx);
sf = [];
stimsizeinsigma = [];
tf = [];
laserdur = [];
stimdur = [];
for ii = 1:numel(idx)
    stro = nex2stro(findfile(char(filename(idx(ii),:))));
    tfidx = strcmp(stro.sum.trialFields(1,:),'tf');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff');
    stimonidx  = strcmp(stro.sum.trialFields(1,:),'stimon_t');
    stimoffidx  = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    lasertrials = logical(stro.trial(:,optstim));
    sf = [sf; stro.sum.exptParams.sf];
    stimsizeinsigma = [stimsizeinsigma; stro.sum.exptParams.sigma];
    tf = [tf; unique(stro.trial(:,tfidx))];
    laserdur = [laserdur; mean(stro.trial(lasertrials,laseroff)-stro.trial(lasertrials,laseron))];
    stimdur = [stimdur; mean(stro.trial(:,stimoffidx)-stro.trial(:,stimonidx))];
end
