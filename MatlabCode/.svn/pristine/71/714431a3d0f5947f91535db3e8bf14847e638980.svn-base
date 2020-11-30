function [modelParams, fv] = LMTFfitting(data, comparison)
    disp('LMTF fitting');
    % This line is just re-representing L- and M-cone contrasts in polar coordinates.
    tf = data(:,3);
    LB = [0   0  1  0 -3 log10(1.00001)];
    UB = [500 1 40 10 -1 log10(2)];
    L = data(:,1);
    M = data(:,2);
    Loog = data(:,4);
    
    initparams = [40 .1 5 3 .005 .002]; % An initial guess for the parameter values
    options = optimset('Display','off');
    if comparison
        thetas = linspace(0,2*pi,6);
    else
        thetas = linspace(0,pi/2,6); 
    end
    thetas(end) = []; toterr = []; fpars = [];
    PLOTINTERMEDIATES = 1;
    % Fitting the data many times (each time rotating the data by a small amount)
    for i = 1:length(thetas)
        rotmat = [cos(thetas(i)) -sin(thetas(i)); sin(thetas(i)) cos(thetas(i))];
        % Rotating data clockwise = rotating fit axis counterclockwise.
        LMrot = [L,M]*rotmat;
        
        % --- GDLH trying new, 1-D initial guesses ---
        % not discounting OOG points, but hopefully this doens't matter
        % because it's just for finding decent initial guesses.
        [th,r] = cart2pol(LMrot(:,1),LMrot(:,2));
        L1 = mod(abs(th),pi) < .17; % a ~20 deg wedge centered on 0
        [fpar1,~] = fmincon(@(params) tf_fiterr(params,data(L1,3),1./r(L1)),initparams,[],[],[],[],LB,UB,[],options);
        L2 =  mod(abs(th-pi/2),pi) < .17; % a ~20 deg wedge centered on pi/2
        [fpar2,~] = fmincon(@(params) tf_fiterr(params,data(L2,3),1./r(L2)),initparams,[],[],[],[],LB,UB,[],options);
        % ---------------------------------------        
        % Now the 2-D fit
        [fpar,fv] = fmincon(@(params) tf_fiterr2(params,[LMrot(:,1) LMrot(:,2) tf],Loog),[fpar1 fpar2],[],[],[],[],[LB LB],[UB UB],[],options);
        fpars(i,:) = fpar;
        toterr(i) = fv;
    end
    bestrotidx = find(toterr == min(toterr));
    %quick hack fix, has only happened once before with utu std inputs 333
    %lum (4 files). Seems to be no difference which bestrotidx is chosen.
    if length(bestrotidx) > 1
       bestrotidx = bestrotidx(1); 
    end
    % Lastly, fitting the theta:
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-8,'Display', 'off'); %used to be "display" "final-detailed"
    [fpar,fv] = fmincon(@(params) tf_fiterr2(params,[L M tf],Loog),[fpars(bestrotidx,1:12) thetas(bestrotidx)],[],[],[],[],[LB LB thetas(bestrotidx)-pi/4],[UB UB thetas(bestrotidx)+pi/4],[],options);
    besttheta = fpar(end);
    modelParams = fpar;
end