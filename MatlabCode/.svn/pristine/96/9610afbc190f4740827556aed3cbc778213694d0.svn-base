function [model, fval, FRsurface, FRpredictedvalues] = fitQuadmodeltospikecounts(WTS,SPC,quadmodelparams)
% Fitting a quadratic model to the spike counts
options = optimset('MaxIter',20000,'MaxFunEvals',30000,'Display','off');
initial_guess = [quadmodelparams max(SPC) 1 0.5 0];

A1 = [-1 0 0 0 0 0 0 0 0;...
    0 -1 0 0 0 0 0 0 0;...
    0 0 0 0 0 -1 0 0 0;...
    0 0 0 0 0 1 0 0 0;...
    0 0 0 0 0 0 -1 0 0;...
    0 0 0 0 0 0 1 0 0;...
    0 0 0 0 0 0 0 -1 0;...
    0 0 0 0 0 0 0 1 0;...
    0 0 0 0 0 0 0 0 -1];
b = [-eps -eps -max(SPC)/2 5*max(SPC) -eps 5 -eps 5 -eps];

[model,fval] = fmincon(@(x)lse(x,WTS,SPC),initial_guess,A1,b,[],[],[],[],[],options); % Then fitting the surface using poisson error model
for ii = 1:100
    [model_tmp,fval_tmp] = fmincon(@(x)lse(x,WTS,SPC),[randn(1,5) max(SPC) randi(4) randi(5) randi(10)],A1,b,[],[],[],[],[],options); % Then fitting the surface using poisson error model
    if fval_tmp < fval
        fval = fval_tmp;
        model = model_tmp;
    end
end
[X,Y] = meshgrid(linspace(-2,2,51));
gensignal = model(1)*X.^2 + model(2)*Y.^2 + model(3)*X.*Y + model(4)*X + model(5)*Y;
FRsurface = (model(6) * (max(gensignal,0).^model(7))./((max(gensignal,0).^model(7) + model(8)^model(7)))) + model(9);

gensignal = model(1)*WTS(:,1).^2 + model(2)*WTS(:,2).^2 + model(3)*WTS(:,1).*WTS(:,2) + model(4)*WTS(:,1) + model(5)*WTS(:,2);
FRpredictedvalues = (model(6) * (max(gensignal,0).^model(7))./((max(gensignal,0).^model(7) + model(8)^model(7)))) + model(9);
    
% least square error
    function err = lse(input,wts,spikecounts)
        A = input(1);
        B = input(2);
        C = input(3);
        D = input(4);
        E = input(5);
        Rmax = input(6);
        N = input(7);
        Chalf = input(8);
        bl = input(9);
        generator_signal = A*wts(:,1).^2 + B*wts(:,2).^2 + C*wts(:,1).*wts(:,2) + D*wts(:,1) + E*wts(:,2);
        mu = (Rmax * (max(generator_signal,0).^N)./((max(generator_signal,0).^N) + Chalf^N)) + bl;
        mu(mu==0) = eps;
        logprob = spikecounts.*log(mu) - mu - log(factorial(spikecounts)); % Poisson error model
        err = -1*sum(logprob);
    end
end

