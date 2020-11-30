function [model,val] = conicsectionfit(x,y,initial_guess)
options.MaxIter = 10000000;
    options.MaxFunEvals = 1000000;
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    
    fval = 10000;
%     [model,val] = fminsearch(@lse, initial_guess, options); % first fitting using LSE
    [model,val] = lsqnonlin(@lse, initial_guess,[],[],options); % using fminunc instead of fminsearch
    
    function err = lse(input)
        A = input(1); 
        B = input(2);
        C = input(3);
        D = input(4);
        E = input(5);
      
        p = A*x.^2 + B*y.^2 + C*x.*y + D*x + E*y - 1;
        err = sum(p.^2);
    end
end

