%% This code resamples the df data to determine a range for each fit variable

function [th_fit_rs,err_rs,rs,Bdist,Cbdist] = resample(datafile,tg,num_resamples)

%% Set Variables
data           = datafile;
xax            = -.4:.001:.4;
p_pred_psych   = nan(num_resamples,length(xax));
p_pred_chron   = p_pred_psych;
th_fit_rs      = nan(num_resamples,length(tg(1:4)));
err_rs         = nan(num_resamples,1);
opts           = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);


%% Reorganize raw data

data(datafile(:,2)==180,1)=data(datafile(:,2)==180,1)*-1;
coh    = data(:,1);
realm  = data(:,2);
choice = data(:,3);
rt     = data(:,4);

listcoh = unique(coh);

for i=1:length(listcoh)
    df(i).coh     = coh(listcoh(i)==coh);
    df(i).sumt    = length(df(i).coh);
    df(i).correct = realm(listcoh(i)==coh)==choice(listcoh(i)==coh);
    df(i).rts     = rt(listcoh(i)==coh);
    df(i).chooseR = choice(listcoh(i)==coh)==0;
end


%% Resample raw data

for n=1:num_resamples
    for i=1:length(listcoh)
        draw                  = randi(length(df(i).coh),length(df(i).coh),1);
        correct               = df(i).correct(draw);
        rts                   = df(i).rts(draw);
        rs(n).df.coh(i,1)     = df(i).coh(1);
        rs(n).df.sumtrials(i) = df(i).sumt;
        rs(n).df.propright(i) = sum(df(i).chooseR(draw))/rs(n).df.sumtrials(i);
        rs(n).df.serc(i)      = std(df(i).chooseR(draw))/rs(n).df.sumtrials(i);
        rs(n).df.meanrt(i)    = mean(rts(correct))*1000;
        rs(n).df.sert(i)      = std(rts(correct)*1000)/sqrt(length(rts(correct)));
    end
end


%% Fit variables 

for z=1:num_resamples
    
    tg(1) = (max(rs(z).df.propright)-min(rs(z).df.propright))/length(rs(z).df.coh);
    tg(2) = sqrt(max(rs(z).df.meanrt));
    tg(3) = min(rs(z).df.meanrt);
    
    [th_fit_rs(z,:),err_rs(z)]=fminsearch('calcparams',tg(1:4),opts,rs(z).df);
    
    k   = th_fit_rs(z,1);
    A   = th_fit_rs(z,2);
    Tnd = th_fit_rs(z,3);
    B   = th_fit_rs(z,4);
    Cb  = B/(2*k*A);
    
    p_pred_psych(z,:) = 1./(1+exp(-2*k*A*(xax-Cb)));
    p_pred_chron(z,:) = (A./(k.*(xax-Cb))).*tanh(k.*A.*(xax-Cb))+Tnd;
    
    
    % Plot results
    figure(6); clf;
    subplot(2,1,1); hold on; grid on;
    axis([-.4 .4 0 1]);
    title('Psychometric Fit to Resampled Data');
    xlabel('%Coherence/Direction');
    ylabel('Proportion Right Choice');
    plot(rs(z).df.coh,rs(z).df.propright,'bo');
    errorbar(rs(z).df.coh,rs(z).df.propright,rs(z).df.serc,'og')
    plot(xax,p_pred_psych(z,:),'b--')
    
    subplot(2,1,2); hold on; grid on;
    title('Chronometric Fit to Resampled Data')
    xlabel('% Coherence / Direction')
    ylabel('Response Time')
    axis([-.4 .4 (.8*min(rs(1).df.meanrt)) (1.2*max(rs(1).df.meanrt))])
    plot(rs(z).df.coh,rs(z).df.meanrt,'go')
    errorbar(rs(z).df.coh,rs(z).df.meanrt,rs(z).df.sert,'go')
    plot(xax,p_pred_chron(z,:),'b--')
    
end

%% Plot distribution of variables (just B for now)

Bdist  = th_fit_rs(:,4);
Cbdist = Bdist./(2*th_fit_rs(:,1).*th_fit_rs(:,2));

figure(25); hold on; grid on;
title('Distribution of B values')
xlabel('B Value')
ylabel('Frequency')
hist(Bdist);

figure(26); hold on; grid on;
title('Distirbution of Cb Values')
xlabel('Cb Value')
ylabel('Occurances')
hist(Cbdist);




