function [model,fval,success] = V1cellfit2(illuminants,surf1,surf2,fundamentals,mean_cone_exc,mode)

hi = size(illuminants,1);
guesses = ones(3,1);
options.MaxIter = 1000;
options.MaxFunEvals = 1000;

[model,fval,success] = fminsearch(@sser,guesses,options);

    function sse = sser(input)
        RGDO_S1 = input;
        RGDO_S2 = -1*RGDO_S1;
        RGDO_S1 = RGDO_S1./norm(RGDO_S1);
        RGDO_S2 = RGDO_S2./norm(RGDO_S2);
        Subunit1_act_cc = [];
        Subunit2_act_cc = [];
        for ii = 1:hi
            net_spectra1 = illuminants(ii,:).*surf1;
            net_spectra2 = (illuminants(ii,:).*surf2);
            
            L1 = net_spectra1 * fundamentals(:,1);
            M1 = net_spectra1 * fundamentals(:,2);
            S1 = net_spectra1 * fundamentals(:,3);
            L2 = net_spectra2 * fundamentals(:,1);
            M2 = net_spectra2 * fundamentals(:,2);
            S2 = net_spectra2 * fundamentals(:,3);
            
            LMS1 = [L1; M1; S1];
            LMS2 = [L2; M2; S2];
            
            % Converting cone excitations into cone contrasts (Weber's
            % contrast)
            LMS1 = (LMS1 - mean_cone_exc)./mean_cone_exc; % cone contrasts
            LMS2 = (LMS2 - mean_cone_exc)./mean_cone_exc; % cone contrasts
            
            Subunit1_act_cc = [Subunit1_act_cc; LMS1'*RGDO_S1];
            Subunit2_act_cc = [Subunit2_act_cc; LMS2'*RGDO_S2];            
        end

%         tmp  = [Subunit1_act_cc Subunit2_act_cc]';
try
    [~,areatmp] = convhull(Subunit1_act_cc,Subunit2_act_cc);
catch
    keyboard;
end
%         crossprod = tmp*tmp';
%         try
%             [~,eigval] = eig(crossprod);
%             ratioeig = max(diag(eigval))/min(diag(eigval));
%             areatmp = ratioeig;
%         catch
%             keyboard;
%         end
        if mode == 1
            sse = (areatmp) + (1-norm(RGDO_S1));
        elseif mode == 2
            sse = (1/areatmp) + (1-norm(RGDO_S1));
        end
    end


end

