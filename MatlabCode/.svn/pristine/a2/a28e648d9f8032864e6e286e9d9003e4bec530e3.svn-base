function V = regularizeSTA(cum_rgbs, cum_n, maxframes, weighted_STA, cosine_bumps, STA)
% Author - Abhishek De, 9/20 
% Finding a regularized STA

options.MaxIter = 1000;
options.MaxFunEvals = 1000;
options.TolFun = 1e-6;
options.TolX = 1e-6;
% maxIter = 1000;
% tolFun = 1e-6;
% options = optimoptions('fmincon','algorithm','trust-region-reflective','gradobj','on','hessian','off','display','iter','maxiter',maxIter,'maxfunevals',maxIter,'tolfun',tolFun,'tolX',tolFun);


guess = [weighted_STA; randn(size(cosine_bumps,2),1)];
lb = -1*ones(size(guess(:)));
ub = 1*ones(size(guess(:)));

% For STA
projs_STA = [];
for ii=1:size(cum_rgbs,2)
    rgbs_STA = cum_rgbs{ii};
    resp_STA = [];
    for jj = 1:size(rgbs_STA,1)
        resp_STA = [resp_STA; conv(rgbs_STA(jj,:),STA(jj,:),'same')];
    end
    projs_STA = [projs_STA; sum(resp_STA,1)'];
end

[xval, fval, success] = fmincon(@(x)err_reg(x), guess(:),[],[],[],[],lb,ub,[], options);
% xval = fminunc(@(x) err_reg(x), guess(:), options);
V = xval(1:300)*(cosine_bumps*xval(301:end))';

    function cost = err_reg(x)
        spatial_weights = x(1:300);
        temporal_weights = x(301:end);
        spatiotemporal_vec = spatial_weights*(cosine_bumps*temporal_weights)';
        x_vec = spatiotemporal_vec(:);
        
        % For spatiotemporal vec 
%         initargs = {x_vec, 0, size(cell2mat(cum_rgbs),2), [100 3 maxframes]};
%         STPROJmod('init',initargs); 
%         for ii = 1:size(cum_rgbs,2)
%             rgbs = cum_rgbs{ii};
%             n = cum_n{ii};
%             STPROJmod(rgbs(:),n);
%         end
%         out = STPROJmod('return');
%         projs = out{1};
%         clear STPROJmod; clear out;
        
        projs_vec = [];
        for aa=1:size(cum_rgbs,2)
            rgbs = cum_rgbs{aa};
            resp = [];
            for bb = 1:size(rgbs,1)
                resp = [resp; conv(rgbs(bb,:),spatiotemporal_vec(bb,:),'same')];
            end
            projs_vec = [projs_vec; sum(resp,1)']; 
        end
        cost  = norm(projs_STA(:)-projs_vec(:),2) + norm(x_vec,2);
        disp(cost);
    end

end

