%% FlashTutorial
% Created Apr_2011 Angueyra
%% Inputs needed
close all
clear all

Stim.Eccentricity = 30; %in degrees of visual angle
Stim.Area = 0.01; %in degrees of visual angle
load('StimSpectra.mat')
Stim.Wavelength=StimSpectra.wavelength;
Stim.Spectra=StimSpectra.blue;
%clear StimSpectra
Stim.Power = 1e-11; %in W (code assumes that this is the power delivered in Stim.Area)
Stim.samplingRate=1e3;
Stim.PreTime=0.5; %in s
Stim.StmTime=1%10*1e-3; %in s
Stim.TailTime= 0.5%-Stim.StmTime; %in s
Stim.DataTime=Stim.PreTime+Stim.StmTime+Stim.TailTime;
Stim.PrePts=Stim.PreTime*Stim.samplingRate;
Stim.StmPts=Stim.StmTime*Stim.samplingRate;
Stim.TailPts= Stim.TailTime*Stim.samplingRate;
Stim.DataPts=Stim.PrePts+Stim.StmPts+Stim.TailPts;
Stim.Stim=zeros(1,Stim.DataPts);
Stim.Stim(Stim.PrePts:Stim.PrePts+Stim.StmPts)=Stim.Power;
TimeAxis=(1:Stim.DataPts)./Stim.samplingRate;

figure(1)
clf
plot(TimeAxis,Stim.Stim,'.-')
%xlabel('Time (s)')
%ylabel('Power (W)')
axis tight
%% Create Cone Grid
% Small grid for this tutorial
Cones.Centers=ConeGrid(Stim.Area,Stim.Eccentricity); %in mm of retina
Cones.CollectingArea = 0.7 * 1e-6; %in mm^2
Cones.Size = 2e-3; %in mm (cell body diameter)
figure(1);clf;hold on;
% MCones
plot(Cones.Centers.MCones.X,Cones.Centers.MCones.Y,'g.')
hold on
for i=1:numel(Cones.Centers.MCones.X)
    %circle([Cones.Centers.MCones.X(i),Cones.Centers.MCones.Y(i)],sqrt(Cones.CollectingArea/pi),20,'g-');
    circle([Cones.Centers.MCones.X(i),Cones.Centers.MCones.Y(i)],Cones.Size/2,20,'g-');
end
% LCones
plot(Cones.Centers.LCones.X,Cones.Centers.LCones.Y,'r.')
for i=1:numel(Cones.Centers.LCones.X)
    %circle([Cones.Centers.LCones.X(i),Cones.Centers.LCones.Y(i)],sqrt(Cones.CollectingArea/pi),20,'r-');
    circle([Cones.Centers.LCones.X(i),Cones.Centers.LCones.Y(i)],Cones.Size/2,20,'r-');
end
% SCones
plot(Cones.Centers.SCones.X,Cones.Centers.SCones.Y,'b.')
for i=1:numel(Cones.Centers.SCones.X)
    %circle([Cones.Centers.SCones.X(i),Cones.Centers.SCones.Y(i)],sqrt(Cones.CollectingArea/pi),20,'b-');
    circle([Cones.Centers.SCones.X(i),Cones.Centers.SCones.Y(i)],Cones.Size/2,20,'b-');
end
hold off
axis tight
axis square
xlabel('Retina (mm)')
ylabel('Retina (mm)')


%% Calibration and Linear Responses

% Normalized Cone Spectral Sensitivities across spectral stimulus wavelength
SpectrumMCone = GeneratePhotoreceptorSpectrum('MCone',Stim.Wavelength);
SpectrumLCone = GeneratePhotoreceptorSpectrum('LCone',Stim.Wavelength);
SpectrumSCone = GeneratePhotoreceptorSpectrum('SCone',Stim.Wavelength);

% Stimulus
RetinalStimArea = VisAngle2Retina_Area(Stim.Eccentricity,Stim.Area);% in mm^2

% Converting from Quantal Power Spectrum to Energy Spectrum
StimEnergySpectra = Stim.Spectra ./ Stim.Wavelength;
% Normalizing total energy to 1
StimEnergySpectra = StimEnergySpectra / trapz(Stim.Wavelength,StimEnergySpectra);
% Making total energy equal to energy measured with power meter (in Watts)
StimEnergySpectra = StimEnergySpectra .* Stim.Power;
% Transforming to photons across wavelengths: Photons(lambda)=Energy(lambda)*lambda/hc
ScFact = 6.6e-34 * 3e8; %Planck's constant * Speed of Light
PhotonSpectra = StimEnergySpectra .* Stim.Wavelength * 1e-9 / ScFact; % (wavelength is in nm)
% Integrating across wavelength to get total number of photons
FluxMCone = sum(PhotonSpectra .* SpectrumMCone);
FluxLCone = sum(PhotonSpectra .* SpectrumLCone);
FluxSCone = sum(PhotonSpectra .* SpectrumSCone);

% Calculating Quantal Catches (Areas need to be in same units)
RStarMCone=  FluxMCone * Cones.CollectingArea / RetinalStimArea; %R*/MCone/s
RStarLCone=  FluxLCone * Cones.CollectingArea / RetinalStimArea; %R*/LCone/s
RStarSCone=  FluxSCone * Cones.CollectingArea / RetinalStimArea; %R*/SCone/s

% Convert Stimulus from W to R*/Cone/s
Stim.StimMConeRStar=(Stim.Stim./max(Stim.Stim)).*RStarMCone;
Stim.StimLConeRStar=(Stim.Stim./max(Stim.Stim)).*RStarLCone;
Stim.StimSConeRStar=(Stim.Stim./max(Stim.Stim)).*RStarSCone;

% Calculate Linear Responses for each Cone Type
LinearResponseMCone=ConeVolution(Stim.StimMConeRStar,Stim.samplingRate);
LinearResponseLCone=ConeVolution(Stim.StimLConeRStar,Stim.samplingRate);
LinearResponseSCone=ConeVolution(Stim.StimSConeRStar,Stim.samplingRate);

figure(4)
subplot (2,2,1)
plot(Stim.Wavelength,Stim.Spectra,'b-')
axis tight
xlabel('Wavelength (nm)')
ylabel('Normalized Stim Power')
title('Stimulus Spectrum')
subplot (2,2,3)
plot(Stim.Wavelength,SpectrumMCone,'g-')
hold on
plot(Stim.Wavelength,SpectrumLCone,'r-')
plot(Stim.Wavelength,SpectrumSCone,'b-')
hold off
axis tight
xlabel('Wavelength (nm)')
ylabel('Normalized Sensitvity')
title('Primate Cone Spectral Sensitivities')
subplot (1,2,2)
plot(TimeAxis,LinearResponseMCone,'g')
hold on
plot(TimeAxis,LinearResponseLCone,'r')
plot(TimeAxis,LinearResponseSCone,'b')
hold off
axis tight
xlabel('Time (s)')
ylabel('Photocurrent (pA)')
title('Linear Responses to 10ms Flash')


%% Adding noise and storing separately responses for each cone type

clear Responses

for i=1:numel(Cones.Centers.MCones.X)
    Noise = ConeNoise(Stim.DataTime,Stim.samplingRate);
    Responses.MCones(i,:)=LinearResponseMCone + Noise;
end

for i=1:numel(Cones.Centers.LCones.X)
    Noise = ConeNoise(Stim.DataTime,Stim.samplingRate);
    Responses.LCones(i,:)=LinearResponseLCone + Noise;
end

for i=1:numel(Cones.Centers.SCones.X)
    Noise = ConeNoise(Stim.DataTime,Stim.samplingRate);
    Responses.SCones(i,:)=LinearResponseSCone + Noise;
end


figure(3)
subplot(4,1,2)
%subplot(3,1,3);hold on;
%plot(TimeAxis,LinearResponseLCone)
%plot(TimeAxis,Responses.LCones(2,:))
plot(TimeAxis,Responses.MCones(2,:),'Color',[0 0.8 0])
hold on
plot(TimeAxis,mean(Responses.MCones),'-','Color',[0 0.5 0])
plot(TimeAxis,LinearResponseMCone,'-.','Color',[0 0.5 0])
hold off
title ('MCone (Example and Mean)')
ylabel('Photocurrent (pA)')

subplot(4,1,1)
plot(TimeAxis,Responses.LCones(1,:),'r')
hold on
plot(TimeAxis,mean(Responses.LCones),'-','Color',[0.5 0 0])
plot(TimeAxis,LinearResponseLCone,'-.','Color',[0.5 0 0])
hold off
title ('LCone (Example and Mean)')
ylabel('Photocurrent (pA)')

subplot(4,1,3)
plot(TimeAxis,Responses.SCones(1,:),'b')
hold on
plot(TimeAxis,mean(Responses.SCones),'-','Color',[0 0.5 0.5])
plot(TimeAxis,LinearResponseSCone,'-.','Color',[0 0.5 0.5])
hold off
title ('SCone (Example and Mean)')
ylabel('Photocurrent (pA)')

subplot(4,1,4)
plot(TimeAxis,Stim.Stim/RetinalStimArea,'k.-')
xlabel('Time (s)')
ylabel('Power (W/mm^2)')