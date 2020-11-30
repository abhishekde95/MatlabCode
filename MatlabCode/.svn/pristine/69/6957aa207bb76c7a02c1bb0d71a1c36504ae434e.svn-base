function Noise=ConeNoise(NoiseLength,samplingRate)
% function Noise=ConeNoise(NoiseLength,samplingRate)
% NoiseLength in s, sampling Rate in Hz
% Function will create cone dark noise based on spectrum
% calculated from recording of dark noise fit with a
% sum of 2 lorentzians up to 1kHz
% Created Jan/2010 Angueyra 
% Values from Jan19_10 (For NIH Grant)
% LorentzCoeffs=[0.0795   60.0    4.0    0.0076  285.0    2.5327];
% Revised Sep/2010 Angueyra
% Using 011310Confocal_c7 and 011310Confocal_c8
% LorentzCoeffs=[0.08   100.0    8.0    0.0076  285.0    2.5327];
% Revised Apr/2011 Angueyra: Changed LorentzCoeffs
% Using 022211Ec07
% Check ConeNoisePSFittingWalkthrough
% LorentzCoeffs=[0.2   30    2.0    0.05  180    2.5];

LorentzCoeffs=[0.2   30    2.0    0.05  180    2.5];

% Number of points to be modelled = NoiseLength*samplingRate
% deltaFreq = 1/NoiseLength
FreqAxis_PS=((1:NoiseLength*samplingRate)-1)./(NoiseLength);
Nyquist=samplingRate/2;
FreqAxis_PS=FreqAxis_PS(1:find(FreqAxis_PS<=Nyquist,1,'last'));
ModelNoisePS=lorentzsum_poles(LorentzCoeffs,FreqAxis_PS);

% PS calculation assumes symmetry of FFT for positive and negative
% frequencies and then squares the magnitude
ModelNoiseFFT_temp = sqrt(ModelNoisePS/2);
% In Matlab's FFT:
%   The first point is the zero-frequency point (DC component), which is also the first point in the PowerSpectrum
%   The following points list the numbers for positive frequencies until Nyquist
%   The second half are the negative frequencies starting at Nyquist
ModelNoiseFFT = [ModelNoiseFFT_temp(1:end) fliplr(ModelNoiseFFT_temp(2:end-1))];
FreqAxisFFT = [-FreqAxis_PS fliplr(FreqAxis_PS(2:end-1))];

% Creating FFT of Gaussian White Noise
WhiteNoiseFFT = fft(normrnd(0,1,1,length(FreqAxisFFT)))*1/samplingRate;


% Shaping White Noise by PowerSpectrum of Cone Dark Noise
ShapedNoiseFFT = complex(ModelNoiseFFT .* real(WhiteNoiseFFT), ModelNoiseFFT .* imag(WhiteNoiseFFT));

ShapedNoise=ifft(ShapedNoiseFFT)*samplingRate;
% Normalize Noise to match original variance
ScaleFactor=sqrt(trapz(FreqAxis_PS,ModelNoisePS)/var(ShapedNoise));
ShapedNoise=ShapedNoise*ScaleFactor;

Noise=real(ShapedNoise);