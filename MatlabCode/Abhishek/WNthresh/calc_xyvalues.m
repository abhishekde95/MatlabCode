function [x_quad,y_quad,rho3] = calc_xyvalues(allthetas, final_model3)
A = final_model3(1); B = final_model3(2); C = final_model3(3); D = final_model3(4); E = final_model3(5); 
rho3 = [];
p = [A*cos(allthetas').^2+B*sin(allthetas').^2+C*(cos(allthetas').*sin(allthetas')) D*cos(allthetas')+E*sin(allthetas') -1*ones(numel(allthetas'),1)];
for kk = 1:size(p,1)
    rts = roots(p(kk,:));
    if all(rts>0)
        r = min(rts); % if both the roots are positive
    else
        r = max(rts);
    end
    rho3 = [rho3; r];
end
L = rho3>0 & rho3==real(rho3);
[x_quad,y_quad] = pol2cart(allthetas(L),rho3(L)');
end

