function hyp = minimize_hyps(hyp, x, y)
if exist(['lbfgsb.' mexext], 'file')
    m = @minimize_lbfgsb;
else
    m = @minimize;
end
try
    hyp = m(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, y);
catch
    m = @minimize;
    hyp = m(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, y);
end
