function [Neurofilter, fval, success] = fitGLM_AD(Stim,response,filterlength)

 % Some options that will be used for function minimization
 error_val = 10000;
 options.MaxIter = 1000000;
 options.MaxFunEvals = 1000000;
 options.TolFun = 1e-6;
 options.TolX = 1e-6;
 
 guesses.stimfilter = randn([filterlength 1]); % initializing the filter
 guesses.postspikefilter = randn([filterlength 1]); % initialzing the post-spike filter
 
 
 [model, fval, success] = fmincon(@minimize_error, [guesses.stimfilter guesses.postspikefilter],[],[],[],[],[],[],[],options); % using gradient-descent based fmincon
 Neurofilter.stimfilter = model(:,1); % Stimulus related filter
 Neurofilter.postspikefilter = model(:,2); % Post-spike filter
 
 % likelihood error function
    function lik = minimize_error(input)
        % MLE assuming Poisson spiking and exponential link function
        filterdrive = conv(Stim,input(:,1),'same'); % Stimulus driven response
        postspikedrive = conv(response,input(:,2),'same'); % History dependent filter
        filterresponse = exp(filterdrive + postspikedrive);
        lik = response.*log(filterresponse) - filterresponse;
        lik = -sum(lik); % maximize log-likelihood ==> minimize negative log-likelihood
    end
end

