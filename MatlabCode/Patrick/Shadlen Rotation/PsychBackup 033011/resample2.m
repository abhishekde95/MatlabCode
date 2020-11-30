%% This code resamples the df data to determine a range for each fit variable

function [th_fit_rs,err_rs] = resample2(datafile,tg,num_resamples)

data=datafile;

data(datafile(:,2)==180,1)=data(datafile(:,2)==180,1)*-1;
coh    = data(:,1);
realm  = data(:,2);
choice = data(:,3);
rt     = data(:,4);

listcoh = unique(coh);

for i=1:length(listcoh)
    df(i).coh     = coh(listcoh(i)==coh);
    df(i).rt      = rt(listcoh(i)==coh);
    df(i).correct = realm(listcoh(i)==coh)==choice(listcoh(i)==coh);
end

for n=1:num_resamples
    for i=1:length(listcoh)
        rs(n).df(i).draw    = randi(length(df(i).coh),length(df(i).coh),1);
        rs(n).df(i).coh     = df(i).coh;
        rs(n).df(i).rt      = df(i).rt(rs(n).df(i).draw);
        rs(n).df(i).correct = df(i).correct(rs(n).df(i).draw);
    end
end


