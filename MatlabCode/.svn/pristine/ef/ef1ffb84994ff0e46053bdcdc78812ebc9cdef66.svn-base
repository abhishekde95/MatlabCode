function monkFits = psychFun(norms, nTrialsByCntrst, a, b, nFigRows, nFigCols, monkFits1,s)

%fits a psychometric function to the choice data and coughs up
%parameter estimates. Works on a case by case basis where each
%"case" is a unique sfs/color combination

%fit the psychometric function
monkFits.performance = monkFits1(1,:);
monkFits.nTrials = monkFits1(2,:);
zeroInd = norms == 0; %don't consider the zero contrast condition
tmpNorms = norms(~zeroInd);
errs = abs(0.82-monkFits.performance(~zeroInd));
aGuess = tmpNorms(find(errs == min(errs), 1, 'last'));
[aSSE, bSSE, ~, success(1)] = weibullFit(norms, monkFits.performance, 'sse', [aGuess 1]);
correctByContrast = (monkFits.performance.*nTrialsByCntrst);
wrongByContrast = (nTrialsByCntrst - correctByContrast);
[monkFits.alpha(a,b), monkFits.beta(a,b), monkFits.gamma(a,b), success(2), modErrs] = weibullFit(norms, [correctByContrast(:), wrongByContrast(:)], 'mle', [aSSE, bSSE]);
dNorm = min(diff(norms))./100;
modelMLE = monkFits.gamma(a,b) + (0.5 - monkFits.gamma(a,b)).*exp(-(([0:dNorm:max(norms)]./monkFits.alpha(a,b)).^monkFits.beta(a,b)));

if all(success);
    %determine errors for the parameter estimates if the SSE and MLE
    %searches have been sucessful
    monkFits.err.mle.alpha(a,b) = modErrs(1);
    monkFits.err.mle.beta(a,b) = modErrs(2);
else
    %nans are used to denote unsucessful fits during subsequent
    %analyses
    monkFits.alpha(a,b) = NaN;
    monkFits.beta(a,b) = NaN;
    monkFits.err.boot.alpha(a,b) = nan;
    monkFits.err.boot.beta(a,b) = nan;
end

figure(6)
subplot(nFigRows,nFigCols,s);hold on;grid on;
h1 = semilogx(norms, monkFits.performance, 'k.');
hold on
h2 = semilogx([0:dNorm:max(norms)], modelMLE, 'k');
hold off
hg = hggroup;  % So that legend works
set(get(get(hg,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
set(h1,'Parent',hg);
set(h2,'Parent',hg);
xlabel('Cone Contrast');
ylabel('p(Correct)');
xlim([(norms(2)./3), max(norms)]); %don't display the zero contrast
ylim([min(monkFits.performance)-.02, 1.02]);

end %psychFun