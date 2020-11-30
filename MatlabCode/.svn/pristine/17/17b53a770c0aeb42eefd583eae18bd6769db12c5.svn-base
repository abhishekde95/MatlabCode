function out = getcumdrive(rsubunit1,rsubunit2,mode)

if mode == 1
    out = (max(0,rsubunit1)).^2 + (max(0,rsubunit2)).^2; % nonlinear integration: half rectified quadratic nonlinearity
elseif mode == 2
    out = sqrt(abs(rsubunit1)) + sqrt(abs(rsubunit2)); % nonlinear integration: full rectified quadratic nonlinearity
elseif mode == 3
    out = sqrt((max(0,rsubunit1))) + sqrt((max(0,rsubunit2))); % nonlinear integration: half rectified quadratic nonlinearity
elseif mode == 4
    out = rsubunit1.^2 + rsubunit2.^2; % full rectified quadratic nonlinearity
elseif mode == 5
    out = rsubunit1 + rsubunit2; % linear integration    
elseif mode ==6 
    out = max(0,(rsubunit1 + rsubunit2).^2); % essentially linear intergration: squaring after adding

end

