uiopen('/Users/greghorwitz/Desktop/MatlabCode/Greg/Sandbox/DataFromSharona.xlsx',1)
% order of columns: 1) WT, 2) WT+Phosp, 3) Mut, 4) Mut+Phosp
ns = sum(~isnan(DataFromSharona));
labelidx = [];
X = [];
for i = 1:length(ns)
    labelidx = [labelidx; mod(i+1,2)*ones(ns(i),1) (i>2)*ones(ns(i),1)];
    X = [X;DataFromSharona(1:ns(i),i)];
end

% Y(:,1) = 0 if no phosp, 1 if phosp
% Y(:,2) = 0 if WT, 1 if Mut


varnames = {'Phosp';'Genotype'};

[P,T,STATS,TERMS] = anovan(log10(X),{labelidx(:,1) labelidx(:,2)},'model','full','varnames',varnames);
multcompare(STATS,'display','on')
%figure; subplot(2,1,1);
%boxplot(DataFromSharona);
%subplot(2,1,2);
%boxplot(log(DataFromSharona));

figure; 
boxplot(log(DataFromSharona));
xlabel('group');
ylabel('log(normalized I)');

subdata = DataFromSharona(1:5,:); % No need to transform data for Friedman test