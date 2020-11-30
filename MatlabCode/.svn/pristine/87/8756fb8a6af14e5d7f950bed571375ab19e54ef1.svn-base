%% ConeDark Noise Power Spectrum Fitting Procedure
% Example code used to derive fit to cone dark noise 
load('StepsExample.mat')
%% Plot Data at Different Light Levels
figure(1)
c=colormap('jet');
TimeAxis=[1:length(StepsExample.Data{1}(1,:))].*StepsExample.samplingInterval;
for i=1:length(StepsExample.Data)
    plot(TimeAxis,StepsExample.Data{i}(1,:),'Color',c(round(i*length(c)/length(StepsExample.Data)),:))
    hold on
    lightLevel{i}=num2str(round(StepsExample.lightLevel{i}));
end
hold off
xlabel('Time(s)')
ylabel('i (pA)')
leg=legend(lightLevel);
set(leg,'box','off','location','best');
%% Plot Power Spectrum for each light level
figure(2)
c=colormap('jet');
loglog(StepsExample.Freq,StepsExample.PS_Dark{i},'.-','Color','k')
hold on
for i=1:length(StepsExample.Data)
    loglog(StepsExample.Freq,StepsExample.PS_Step{i},'.-','Color',c(round(i*length(c)/length(StepsExample.Data)),:))
end
hold off
xlabel('Freq(Hz)')
ylabel('Power (pA^2/Hz)')
leg=legend(lightLevel);
set(leg,'box','off','location','best');
xlim([1 5000])
ylim([1e-5 5])
%% Calculating Cone Dark Noise
% Assume that Brightest Step is Instrumental Noise. Real Dark Noise = PS_Dark - PS_BrightestStep
figure(3)
loglog(StepsExample.Freq,StepsExample.DiffPS,'.-')
hold all
loglog(StepsExample.Freq_LinearFilter,StepsExample.PS_LinearFilter*5e4,'k-')
hold off
xlabel('Freq(Hz)')
ylabel('Power (pA^2/Hz)')
xlim([1 5000])
ylim([1e-7 1])
leg=legend('DarkNoisePS','SinglePhotonResponsePS(scaled)');
set(leg,'box', 'off','location','best')
%% Fitting Cone Dark Noise
% Can't get nlinfit to converge properly, but fit by eye seems good enough
LorentzCoeffs=[0.16   55    4    0.045  190    2.5];
Fit=lorentzsum_poles(LorentzCoeffs,StepsExample.Freq);
Beta=nlinfit(StepsExample.DiffPS(StepsExample.DiffPS>0),StepsExample.Freq(StepsExample.DiffPS>0),@lorentzsum_poles,LorentzCoeffs);
disp(Beta)
BetaFit=lorentzsum_poles(Beta,StepsExample.Freq);
figure(4)
loglog(StepsExample.Freq(StepsExample.DiffPS>0),StepsExample.DiffPS(StepsExample.DiffPS>0),'b.-')
hold all
loglog(StepsExample.Freq_LinearFilter,StepsExample.PS_LinearFilter*0.5e5,'k-')
loglog(StepsExample.Freq,Fit,'r-','LineWidth',1.5)
loglog(StepsExample.Freq,BetaFit,'go-')
hold off
xlim([1 1000])
ylim([5e-4 0.5])
xlabel('Freq(Hz)')
ylabel('Power (pA^2/Hz)')