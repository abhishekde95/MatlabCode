function [LFPmean, LFPsem, BASEmean, BASEsem, RESPmean, RESPsem] = SMurray_freq_spectrum(near, far, nearBase, farBase, params)
                                                                           
                                                                           
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
            [sLFPnear,fLFPnear] = mtspectrumc(LFPnear, params);
            [sLFPfar,fLFPfar] =   mtspectrumc(LFPfar, params);
            
            %LFP power during fixation, prior to stimulus presentation
            BASEnear = nearBase(k,:,j,i);
            BASEfar  =  farBase(k,:,j,i);
            [sBASEnear,fBASEnear] = mtspectrumc(BASEnear, params);
            [sBASEfar,fBASEfar] =   mtspectrumc(BASEfar, params);
            
            %{
            %plot for range of frequencies in 'params.fpass'
            figure; 
            plot(fLFPnear, 10*log10(sLFPnear))
            set(gca,'FontName','Times New Roman','Fontsize', 16);
            xlabel('Frequency (Hz)'); ylabel('Spectrum (dB)');
            title('LFP spectrum');
            %}
            
            %convert to dB
            tmpnearL(k,1) = 10*log10(sLFPnear); 
            tmpfarL(k,1) = 10*log10(sLFPfar); 
            tmpnearB(k,1) = 10*log10(sBASEnear);
            tmpfarB(k,1) = 10*log10(sBASEfar); 
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
