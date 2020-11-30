function LinearResponse=ConeVolution(Stim,samplingRate)
% function LinearResponse=ConeVolution(Stim,samplingRate)
% Convolution between Stimulus (in R*/s) and Cone Linear Filter (in pA/R*)
% LinearFilter derived from fitting linear filter derived from white noise
% stimulation of L and M cones (check ConeFilterFittingWalkthrough)
% Created Apr_2011 Angueyra

TimeAxis=1:length(Stim);
TimeAxis=TimeAxis./samplingRate;
FilterCoeffs=[0.6745    0.0216    0.0299    0.5311   34.1814];
Filter=ConeEmpiricalDimFlash(FilterCoeffs,TimeAxis);

FilterFFT=fft(Filter)./samplingRate;
StimFFT=fft(Stim)./samplingRate;
LinearResponse=ifft(FilterFFT .* StimFFT) .* samplingRate;

% figure(2)
% subplot(3,1,1)
% plot(TimeAxis,Stim,'.-')
% subplot(3,1,2)
% plot(TimeAxis,Filter,'.-')
% 
% %xlim([0 0.3])
% subplot(3,1,3)
% plot(TimeAxis,LinearResponse,'.-')
