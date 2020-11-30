function [AUC] = HorwitzRiekeROC(Responses)

% Variables
sig = nan(3,3,size(Responses.LCones,1));
AUC = nan(size(Responses.LCones,1),1);

mu  = [mean(Responses.LCones,2),mean(Responses.MCones,2),mean(Responses.SCones,2)];

% Calcuate ROC
for s = 1:size(Responses.LCones,1) % March through all 15 stimuli
    
    sig(:,:,s) = cov([Responses.LCones(s,:);Responses.MCones(s,:);Responses.SCones(s,:)]');
    
    if s==1
        w = (sig(:,:,1)^-1)*mu(1,:)';
    else
        w = ((sig(:,:,1) + sig(:,:,s))^-1) * (mu(s,:)' - mu(1,:)');
    end
    
    projection0 = sort([Responses.LCones(1,:); Responses.MCones(1,:); Responses.SCones(1,:)]' * w);
    projection1 = sort([Responses.LCones(s,:); Responses.MCones(s,:); Responses.SCones(s,:)]' * w);
    
    thresh = unique([projection0 projection1]);
    
    TPR = nan(1,length(thresh));
    FPR = TPR;
    
    for n=1:length(thresh)
        TPR(n) = sum(projection0 >= thresh(n))/length(projection0);
        FPR(n) = sum(projection1 >= thresh(n))/length(projection1);
    end
    
    
    figure(2);hold on; grid on;
    if s==1
        clf;
    end
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    axis([-.1 1.1 -.1 1.1])
    h = plot(TPR,FPR);
    set(h,'Color',unifrnd(0,1,1,3));
    
    AUC(s) = trapz(fliplr(TPR),fliplr(FPR));
    
%     figure(4);clf;hold on;grid on;
%     plot3(Responses.LCones(1,:),Responses.MCones(1,:),Responses.SCones(1,:),'bo')
%     plot3(Responses.LCones(s,:),Responses.MCones(s,:),Responses.SCones(s,:),'go')
%     plot3(mu(1,1),mu(1,2),mu(1,3),'k*')
%     plot3(mu(s,1),mu(s,2),mu(s,3),'k*')
%     plot3([0;w(1)*1e14],[0;w(2)*1e14],[0;w(3)*1e14],'k--')
%     plot3([mu(1,1); mu(1,1)+1e13*w(1)],[mu(1,2); mu(1,2)+1e13*w(2)],[mu(1,3); mu(1,3)+1e13*w(3)],'k-')
%     
end

