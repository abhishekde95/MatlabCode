function [LFPmean, LFPsem, BASEmean, BASEsem, RESPmean, RESPsem] = SMurray_freq_filter(near, far, ...
                                                                                       nearBase, farBase, filtkernel1, filtkernel2)

%analyze each stim type, channel, and trial separately (SEM based on trials)
LFPmean = NaN(2, size(near,4), size(near,3));
LFPsem =  NaN(2, size(near,4), size(near,3));
BASEmean = NaN(2, size(near,4), size(near,3));
BASEsem =  NaN(2, size(near,4), size(near,3));
RESPmean = NaN(2, size(near,4), size(near,3));
RESPsem =  NaN(2, size(near,4), size(near,3));

for i = 1:size(near,4) %stimulus radius
    for j = 1:size(near,3) %channel
        tmpnearL = NaN(size(near,1), 1);
        tmpfarL  = NaN(size(near,1), 1);
        tmpnearB = NaN(size(near,1), 1);
        tmpfarB  = NaN(size(near,1), 1);
        
        for k = 1:size(near,1) %trial
            
            %LFP power during stimulus presentation
            LFPnear = near(k,:,j,i);
            LFPfar  =  far(k,:,j,i);
            filteredLFPnear = sum(conv(LFPnear, filtkernel1,'same').^2+conv(LFPnear, filtkernel2,'same').^2);
            filteredLFPfar = sum(conv(LFPfar, filtkernel1,'same').^2+conv(LFPfar, filtkernel2,'same').^2);

            %LFP power during fixation, prior to stimulus presentation
            BASEnear = nearBase(k,:,j,i);
            BASEfar  =  farBase(k,:,j,i);
            filteredBASEnear = sum(conv(BASEnear, filtkernel1,'same').^2+conv(BASEnear, filtkernel2,'same').^2);
            filteredBASEfar = sum(conv(BASEfar, filtkernel1,'same').^2+conv(BASEfar, filtkernel2,'same').^2);
           
            tmpnearL(k,1) = filteredLFPnear;
            tmpfarL(k,1) = filteredLFPfar;
            tmpnearB(k,1) = filteredBASEnear;
            tmpfarB(k,1) = filteredBASEfar;
        end
        
        %calculate LFP response (LFP power during stim - LFP power prior to stim)
        tmpnearR = tmpnearL - tmpnearB; 
        tmpfarR = tmpfarL - tmpfarB;
        
        %Mean and SEM for LFP power during stimulus, across trials
        LFPmean(1,i,j) = mean(tmpnearL);
        LFPmean(2,i,j) = mean(tmpfarL);
        LFPsem(1,i,j) = std(tmpnearL) / sqrt(length(tmpnearL));
        LFPsem(2,i,j) = std(tmpfarL) / sqrt(length(tmpfarL));
        
        %Mean and SEM for LFP power prior to stimulus (baseline), across trials
        BASEmean(1,i,j) = mean(tmpnearB);
        BASEmean(2,i,j) = mean(tmpfarB);
        BASEsem(1,i,j) = std(tmpnearB) / sqrt(length(tmpnearB));
        BASEsem(2,i,j) = std(tmpfarB) / sqrt(length(tmpfarB));
        
        %Mean and SEM for LFP response (change-from-baseline), across trials
        RESPmean(1,i,j) = mean(tmpnearR);
        RESPmean(2,i,j) = mean(tmpfarR);
        RESPsem(1,i,j) = std(tmpnearR) / sqrt(length(tmpnearR));
        RESPsem(2,i,j) = std(tmpfarR) / sqrt(length(tmpfarR));
        
    end
end


end