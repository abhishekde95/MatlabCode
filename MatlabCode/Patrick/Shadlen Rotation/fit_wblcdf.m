function [final_params,final_error, eflag] = fit_wblcdf(x,y)

% setting a bunch of fitting options
fittingOptions = optimset('tolfun', .00001, ...
                            'tolx', .01, ...
                            'diffmaxchange', .001, ...
                            'display', 'off', ...
                            'maxfunevals', 4000, ...
                            'maxiter', 3000, ...
                            'RelLineSrchBnd', 1e-1, ...
                            'RelLineSrchBndDuration', 20, ...
                            'largeScale', 'off');
fittingOptions = optimset(@fminsearch);




% initial guesses for alpha/beta/gamma
% alpha is the position of threshold
% beta is the slope
% gamma is the shift of the function (im not guessing it right now)
init= [2,1];

[final_params,final_error, eflag] = fminsearch(@(params) fit_error(x,y,params),...
    init, fittingOptions);



    function ferror = fit_error(x_data,y_data,params)
        alpha   = params(1);        % position
        beta    = params(2);        % slope
      %  gamma   = params(3);        % shift

        if (alpha<=0)| (beta<=0)
            ferror = 1e10;
            return
        end
%         if (gamma>0)
%             error = 1e10;
%             return;
%         end

        % if you want to include gamma use this funtion
        %expected_p = 1- exp(-((x_data - gamma)/alpha).^beta);

        
        expected_p=  1 - 0.5*(exp(-(x_data/alpha).^beta));
        
        
        % sanity checks - make error very large is the values 
        % returned are non-sensical
        expected_p(expected_p <1e-4) = 1e-12;
        expected_p(expected_p>(1-1e-4)) = 1 - 1e-12;

        % actual data
        observed_p = y_data;
        observed_p(observed_p <1e-4) = 1e-12;
        observed_p(observed_p>(1-1e-4)) = 1 - 1e-12;

        expected_q = 1- expected_p;
        observed_q = 1- observed_p;


        observed_mle= 0;
        expected_mle= 0;
        for i=1:length(x_data)

            observed_mle = observed_mle + 1000*log((observed_p(i)^observed_p(i))*(observed_q(i)^observed_q(i)));
            expected_mle = expected_mle + 1000*log((expected_p(i)^observed_p(i))*(expected_q(i)^observed_q(i)));

        end
        
        clf;
        semilogx(x_data, y_data, 'k.','MarkerSize', 25); hold on;
        semilogx(x_data, expected_p,'r'); 
        set(gca, 'ytick',0.5:0.1:1)
        axis([0 100 0.5 1.05])
        xlabel('coherence(%)');
        ylabel('proprtion correct');
        pause(.5);
        
        
        ferror = 2*(observed_mle - expected_mle);

    end
end
