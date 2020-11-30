% This code is meant to return voltage at a single point through time
% JPW May_2011
% This code is rendered obsolete by stimDynamics.m

function [stimvals] = gabordynamics(nstd,sigma,nsindur,nsinundgaus,t_total)
% nstd = # of standard deviations in gaussian
% sigma = pixels per std
% nsindur = # of sinusoid cycles for duration of stimulus
% nsinundgaus = # of sinusoid cycles under the gaussian at one time

% NOTE: nstd*sigma should equal nsin*sigsin
% NOTE: units of space should be pixels, time should be ms

% Spatial Parameters
npixels = nstd*sigma; %Diameter of stim in pixels
x_rad  = linspace(0,nsinundgaus*2*pi,npixels); %Pixels converted to radians
x = linspace(-nstd/2*sigma, nstd/2*sigma, npixels); %Spatial coordinates for gaussian

% Temporal Parameters
dx = (nsindur*2*pi)/t_total;
tpts = 0:t_total;

gaussian = normpdf(x,0,sigma);
gaussian = gaussian./max(gaussian);

stimvals = nan(numel(gaussian),numel(gaussian),numel(tpts));
disp(numel(tpts));
for t = 1:numel(tpts)
    t
    for s = 1:numel(gaussian)
        sinusoid = sin(x_rad+(dx*tpts(t)));
        stimvals(s,:,t) = sinusoid.*gaussian.*gaussian(s);
        %figure(1);clf;hold on;
        %plot(x,gaussian.*sinusoid)
        
%         subplot(2,1,1)
%         cla,hold on,
%         plot(x, gaussian, 'k');
%         plot(x, sinusoid, 'r');
%         hold off
%         
%         subplot(2,1,2)
%         plot(x, gaussian .* sinusoid)
%         ylim([-1 1])
    end
end