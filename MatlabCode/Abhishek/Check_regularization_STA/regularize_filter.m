function V = regularize_filter(stim, output, guess, cosine_bumps, mode)

% mode == 1, normal regularization
% mode == 2, regularization with L1 minimization
% mode == 3, regularization with L2 minimization

% A function for regularizing the filter 
options.MaxIter = 1e4;
options.MaxFunEvals = 1e4;
options.TolFun = 1e-6;
options.TolX = 1e-6;

lb = -1*ones(size(guess(:)));
ub = 1*ones(size(guess(:)));
lambda = 1; 

if mode == 1
    % Normal regularization
    [V, fval, success] = fmincon(@(x)err_reg(x), guess(:),[],[],[],[],lb,ub,[], options);
    
elseif mode == 2
    % Regularization with L1 norm
    [V, fval, success] = fmincon(@(x)err_reg_L1(x), guess(:),[],[],[],[],lb,ub,[], options);

    
elseif mode == 3
    % Regularization with L2 norm
    [V, fval, success] = fmincon(@(x)err_reg_L2(x), guess(:),[],[],[],[],lb,ub,[], options);

elseif mode ==4 
    % Regularization using Greg's method 
    STA = guess;
    [xval, fval, success] = fmincon(@(x)err_reg_GH(x), randn(size(cosine_bumps,2),1),[],[],[],[],lb,ub,[], options);
    
    if size(xval,2)~=1
        xval = xval';
    end
    V = cosine_bumps*xval;
    
end

    function cost = err_reg(x)
        % Normal regularization   
        projs = conv(stim,x,'same'); % Projecting stimulus onto the recovered filter
        cost = min([rocN(projs(output==0),projs(output>0)) 1-rocN(projs(output==0),projs(output>0))]); % similar to MID approach 

    end

    function cost = err_reg_L1(x)
        % Normal regularization   
        projs = conv(stim,x,'same'); % Projecting stimulus onto the recovered filter
        cost = min([rocN(projs(output==0),projs(output>0)) 1-rocN(projs(output==0),projs(output>0))]); % similar to MID approach 

    end

    function cost = err_reg_L2(x)
        % Regularization with L2 norm 
        projs = conv(stim,x,'same'); % Projecting stimulus onto the recovered filter
        cost = min([rocN(projs(output==0),projs(output>0)) 1-rocN(projs(output==0),projs(output>0))]) + lambda*norm(x(:),2); % similar to MID approach

    end

    function cost = err_reg_GH(x)
        % Regularization with L2 norm using Greg's method 
        vec = cosine_bumps*x;
        projs = conv(stim,vec','same'); % Projecting stimulus onto the recovered filter
        projs_STA = conv(stim,STA,'same');
        cost = norm(projs(:)-projs_STA(:),2) + lambda*norm(vec(:),2); % similar to MID approach

    end
end

