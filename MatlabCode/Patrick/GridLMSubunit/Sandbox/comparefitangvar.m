global GLMSPopData


GLMSPopAnalysis
datanames = GLMSPopData(1,:);
surfparams = [GLMSPopData{2:end,strcmp(datanames,'Surface Parameters')}];
oneDL = strcmp(GLMSPopData(2:end,strcmp(datanames,'nD')),'1D');
twoDL = strcmp(GLMSPopData(2:end,strcmp(datanames,'nD')),'2D');

% Variance on preferred direction
oneDsurfs = [surfparams(oneDL).oneD];
twoDsurfs = [surfparams(twoDL).twoD];
oneDdirvar = nan(sum(oneDL),1);
twoDdirvar = nan(sum(twoDL),1);
for n = 1:numel(oneDdirvar)
    oneDdirvar(n) = sqrt(1/oneDsurfs(n).Hessian(end-1,end-1)/pi*180);
end
for n = 1:numel(twoDdirvar)
    twoDdirvar(n) = sqrt(1/twoDsurfs(n).Hessian(end-1,end-1)/pi*180);
end  

figure(1); clf; hold on; 
hist(oneDdirvar,50)

figure(2); clf; hold on;
hist(twoDdirvar,50,'r')

