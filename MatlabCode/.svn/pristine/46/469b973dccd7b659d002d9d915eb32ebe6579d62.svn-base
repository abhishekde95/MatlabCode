function [df]=orgrawdata(datafile,q)

%% Reorganize and relabel data
% Label leftward choices as negative
% Convert: right choice = 0, left choice = 1
datafile(datafile(:,2)==180,1)=datafile(datafile(:,2)==180,1)*-1;
%datafile(datafile(:,3)==180,3)=1;
coherence=datafile(:,1);
df.coh=unique(coherence);
realmotion=datafile(:,2);
choice=datafile(:,3);
reactiontimes=datafile(:,4);


%% Separate datafiles by columns
% Then organize by coherence
% Get total trials at each coh, proportion of rightward choices, correct choices, se, 
% mean reaction times, s
for i=1:length(df.coh)
    L=coherence==df.coh(i);
    df.sumtrials(i)=sum(L);
    df.propright(i)=sum(choice(L)==0)/df.sumtrials(i);
    df.serc(i)=std(choice(L)==0)/sqrt(df.sumtrials(i));
    rts=reactiontimes(L);
    df.meanrt(i)=mean(rts(realmotion(L)==choice(L)))*1000;
    df.stdrt(i)=std(rts(realmotion(L)==choice(L)))*1000;
    df.sert(i)=std(rts(realmotion(L)==choice(L))*1000)/sqrt(length(rts(realmotion(L)==choice(L))));
end 


%% Plot results    
figure(q); clf;
subplot(2,1,1); hold on; grid on;
axis([-.4 .4 0 1])
xlabel('% Coherence/Direction')
ylabel('Proportion Right Choice')
title('Raw Psychometric Data')
errorbar(df.coh,df.propright,df.serc,'og')
plot(df.coh,df.propright,'go')

subplot(2,1,2); hold on; grid on;
axis([-.4 .4 .9*min(df.meanrt) 1.1*max(df.meanrt)])
xlabel('% Coherence/Direction')
ylabel('Reaction Time')
title('Raw Chronometric Data')
errorbar(df.coh,df.meanrt,df.sert,'go')
plot(df.coh,df.meanrt,'go')
