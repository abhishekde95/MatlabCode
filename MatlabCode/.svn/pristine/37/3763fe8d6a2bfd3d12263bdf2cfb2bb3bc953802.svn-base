
function DerivCheck(funptr, X0, opts, varargin)
% DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
%
%  Checks the derivative of a function 'funptr' at a point X0, 
%  for purposes of optimization
%
%  Call with same arguments as you would call the optimization routine.


tol = 1e-8;  % Size of numerical step to take
rr = randn(length(X0),1)*tol;  % Generate small random-direction vector

[~,JJ] = feval(funptr,X0,varargin{:});

f1 = feval(funptr, X0-rr/2, varargin{:});  % Evaluate function at X0 - eps
f2 = feval(funptr, X0+rr/2, varargin{:});  % Evaluate function at X0 + eps

fprintf('Derivs: Analytic vs. Finite Diff = [%.4e, %.4e]\n', dot(rr, JJ), f2-f1);
