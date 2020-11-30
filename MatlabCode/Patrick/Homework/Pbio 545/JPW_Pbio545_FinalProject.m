% This script models the recreation of a V1 receptive field from Fourier
% analysis of impulse response functions of a model V1 cell.

% 3/12  Created.    JPW

clear all
%close all

% In this script, I attempt to recreate a V1 cell from some of its
% component processing areas.  I begin with a stimulus, which is defined in
% terms of energy at particular wavelengths in the visible light spectrum.
% From there, I calculate how many photons are collected by each type of 
% cone photoreceptor, and I follow the signal cascade inside of each cone 
% type to determine the change of current that results.  I use that
% change current as input into a retinal ganglion cell, which has both
% excitatory and suppressive regions (on center, off surround).  Then,
% wiring several adjacent ganglion cells together, I create the receptive
% field for a V1 cell.  This cell should mimick the activity of a simple V1
% cell, and should also be suceptible to the same analysis as a typical V1
% cell.  I will attempt to recreate the receptive field of my model V1 cell
% through the same analysis used by Movshon et al in their 1977 papers.


%% Part I: Stimulus

% Define a stimulus.  We begin by determining the output of each individual
% cone type to a full field flash with specified component frequencies at
% a specified power.  Since I will be ultimately calculating the resulting
% current of a photoreceptor as the cube of cGMP concentration (a non
% linearity!), I want to keep the change in cGMP fairly small - that is,
% small enough that the change is current is roughly linear.  Thus, we want
% to define the power as some small value.

% Here, we can define several parameters of the stimulus.  We must first
% specify the time course of the stimulus.  We will model the cones as
% sitting in the dark, then recieving a dim flash of light.  First, we
% specify how much time we observe the cone before and after the flash, as
% well as the duration of the flash itself.  I am choosing 64ms to mimic
% the experimental conditions of Movshon et al.

%%%% User-Defined Parameters %%%%

% Temporal parameters of the stimulus
Stim.PreTime = .5; % in s
Stim.PostTime = 2; % in s
Stim.StimTime = .64; % in s
Stim.SamplingRate = 1000; % samples/sec

% Spectrum of wavelengths relevant to the cones.  (Must define the stimulus
% over this range of  wavelengths and at 5nm intervals for now.  This is 
% because I'm using cone spectral sensitivites calculated elsewhere and
% importing them from a mat file.)
Stim.Wavelength = 350:5:800; % in nanometers 

% Choose the relative energy spectrum of the stimulus (Note: this is a 
% normalized intensity at each wavelength which will be scaled to integrate 
% to specified power of the flash).
Stim.Spectrum = zeros(size(Stim.Wavelength));
Stim.Spectrum(find(Stim.Wavelength == 400):find(Stim.Wavelength == 500)) = 1;
Stim.Spectrum = Stim.Spectrum ./ sum(Stim.Spectrum);

% Choose total power of the stimulus in a given area (W / mm.^2 on the retina))
Stim.Power = 1e-9; % in J/s
RetinalStimArea = .01; % in mm^2


%%%% Calculations %%%%

% Setting up stimulus temporal parameters
Stim.PrePts = Stim.PreTime * Stim.SamplingRate; % in samples
Stim.StimPts = Stim.StimTime * Stim.SamplingRate;
Stim.PostPts = Stim.PostTime * Stim.SamplingRate;
Stim.DataPts = Stim.PrePts + Stim.StimPts + Stim.PostPts;
TimeAxis = (1:Stim.DataPts) ./ Stim.SamplingRate; % in s

% Square Pulse
Stim.Stim = zeros(1,Stim.DataPts);
% Stim.tStim = sin((2*pi*Stim.NCycles/Stim.StimPts):...
%     (2*pi*Stim.NCycles/Stim.StimPts):2*pi*Stim.NCycles);
% Stim.tStim(sign(Stim.tStim)==1) = 1; %Normalized value (Scaled by R*'s later)
% Stim.tStim(sign(Stim.tStim)==-1) = 0;
% Stim.Stim(Stim.PrePts+1:Stim.PrePts+Stim.StimPts) = Stim.tStim;
Stim.Stim(Stim.PrePts+1:Stim.PrePts+Stim.StimPts) = 1;

% Square Grating
% Stim.Stim = ones(1,Stim.DataPts) * Stim.Power/2;
% Stim.tStim = sin((2*pi*Stim.NCycles/Stim.StimPts):...
%      (2*pi*Stim.NCycles/Stim.StimPts):2*pi*Stim.NCycles);
% Stim.tStim(sign(Stim.tStim)==1) = 1;
% Stim.tStim(sign(Stim.tStim)==-1) = 0;
% Stim.Stim(Stim.PrePts+1:Stim.PrePts+Stim.StimPts) = Stim.tStim;

% Sinusoidal Grating
% Stim.Stim = ones(1,Stim.DataPts) * Stim.Power/2;
% Stim.tStim = ((sin((2*pi*Stim.NCycles/Stim.StimPts):...
%     (2*pi*Stim.NCycles/Stim.StimPts):2*pi*Stim.NCycles))+1)./2;
% Stim.Stim(Stim.PrePts+1:Stim.PrePts+Stim.StimPts) = Stim.tStim;


% Converting from Quantal Power Spectrum to Energy Spectrum (Energy at each wavelength)
StimEnergySpectra = Stim.Spectrum ./ Stim.Wavelength;
% Normalizing total energy to 1
StimEnergySpectra = StimEnergySpectra ./ sum(StimEnergySpectra);
% Making total energy equal to specified stimulus power (in Watts)
StimPowerSpectra = StimEnergySpectra .* Stim.Power;
% Transforming to photons across wavelengths: Photons(lambda)=Energy(lambda)*lambda/hc
ScFact = 6.6e-34 * 3e8; %Planck's constant * Speed of Light
PhotonSpectra = StimPowerSpectra .* Stim.Wavelength * 1e-9 / ScFact; % (wavelength is in nm)


% Plot stimulus
figure(1); clf; hold on;
subplot(3,1,1)
plot(Stim.Wavelength,Stim.Spectrum,'*-');
title('Stimulus Spectrum')
xlabel('Wavelength (nm)')
ylabel('Relative Intensity')
subplot(3,1,2)
plot(Stim.Wavelength,PhotonSpectra,'*-')
title('Photons Delivered by Stimulus')
xlabel('Wavelength (nm)')
ylabel('Photons/s/mm^2')
subplot(3,1,3)
plot(TimeAxis,Stim.Stim,'o-')
axis tight
title('Stimulus Time Course')
xlabel('Time (s)')
ylabel('Watts (J/s)')


%% Part II: Converstion to R*'s

% We now know the number of photons at each wavelength that is being
% emitted from the stimulus.  We will here ignore the complications of
% losing photons due to confounding factors like pre retinal filtering.
% Instead, we will just pretend that we can specify the stimulus parameters
% at some unit area of the retina.  Here, we will scale the stimlus by the
% size of a cone to calculate how many of these photons/s/unit area are
% falling on the photoreceptors.  Since each photoreceptor has a different
% absorbtion spectrum, and because we assigned our stimulus with a narrow 
% band of wavelengths, each of the cones will absorb a different amount of
% photons/s.  If each photons activates a single rhodopsin, photons/s
% should be equal to R*s/s.

% Size constants for the cones
Cones.CollectingArea = 0.7 * 1e-6; %in mm^2
Cones.Size = 2e-3; %in mm (cell body diameter)

% Load Normalized Cone Spectral Sensitivities across spectral stimulus wavelength
load SpectralSensitivities.mat
SpectrumLCone = SpectralSensitivites(:,1);
SpectrumMCone = SpectralSensitivites(:,2);
SpectrumSCone = SpectralSensitivites(:,3);

% Integrating across wavelength to get total number of photons
FluxLCone = sum(PhotonSpectra .* SpectrumLCone');
FluxMCone = sum(PhotonSpectra .* SpectrumMCone');
FluxSCone = sum(PhotonSpectra .* SpectrumSCone');

% Calculating Quantal Catches (Photons falling on  a single cone)
PhotonsLCone =  FluxLCone * Cones.CollectingArea / RetinalStimArea; %R*/LCone/s
PhotonsMCone =  FluxMCone * Cones.CollectingArea / RetinalStimArea; %R*/MCone/s
PhotonsSCone =  FluxSCone * Cones.CollectingArea / RetinalStimArea; %R*/SCone/s
Cones.Photons = cat(1,PhotonsLCone,PhotonsMCone,PhotonsSCone);

figure(2); clf; hold on;
plot(Stim.Wavelength,SpectrumLCone,'r')
plot(Stim.Wavelength,SpectrumMCone,'g')
plot(Stim.Wavelength,SpectrumSCone,'b')
plot(Stim.Wavelength,Stim.Spectrum./max(Stim.Spectrum),'k--')
ylim([-.1 1.1])
title('Spectral Sensitivities and Stimulus Spectrum')
xlabel('Wavelength (nm)')
ylabel('Relative Sensitivity/Intensity')
legend('L Cone','M Cone','S Cone','Stimulus')


%% Part III: Model Single Cone Response

% Here we calculate the change in current in each cell as a result of 
% our flashed stimulus.


% Rate Constants, etc.
TimeStep = 1/Stim.SamplingRate;
NumPts = numel(Stim.Stim);
sigma = 5;				% rhodopsin activity decay rate constant (1/sec)
phi = 5;				% phosphodiesterase activity decay rate constant (1/sec)
eta = 10;				% phosphodiesterase activation rate constant (1/sec)
gdark = 15;				% concentration of cGMP in darkness
cgmp2cur = 8e-3;		% constant relating cGMP to current
cdark = 0.5;			% dark calcium concentration
beta = 2;				% rate constant for calcium removal in 1/sec
hillcoef = 4;			% cooperativity
hillaffinity = 0.3;		% affinity

cur2ca = beta * cdark / (cgmp2cur * gdark^3);				% get q using steady state
smax = eta/phi * gdark * (1 + (cdark / hillaffinity)^hillcoef);		% get smax using steady state

figure(3); clf; hold on;
Cones.Responses = nan(3,numel(Stim.Stim));

% Loop through the three cone types
for n = 1:3
    
    % Tailor stimulus to each cone type
    r = Stim.Stim .* Cones.Photons(n,:)/Stim.SamplingRate;
    
    % Initial Conditions
    ar(1) = 0;
    p(1) = eta/phi;
    g(1) = gdark;							% initial condition
    s(1) = gdark * eta/phi;		% steady state: synthesis = hydrolysis
    c(1) = cdark;
    
    % solve difference equations
    for pnt = 2:numel(Stim.Stim)
        ar(pnt) = ar(pnt-1) + (r(pnt-1) - sigma * ar(pnt-1)) * TimeStep;
        p(pnt) = p(pnt-1) + TimeStep * (ar(pnt-1) + eta - phi * p(pnt-1));
        g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
        c(pnt) = c(pnt-1) + TimeStep * (cur2ca * cgmp2cur * g(pnt-1)^3 - beta * c(pnt-1));
        s(pnt) = smax / (1 + (c(pnt) / hillaffinity)^hillcoef);
    end
    cur = cgmp2cur * g.^3;
    
    Cones.Responses(n,:) = cur;
    
    line = ['r';'g';'b'];
    
    % Plot
    figure(3); hold on; nf = 7; q = 1;
    subplot(nf,1,q); q=q+1; hold on;
    plot(TimeAxis,r,['-' num2str(line(n))])
    %axis tight
    title('Photons to Current')
    ylabel('Absorbed Photons')
    legend('L Cones','M Cones','S Cones','Location','NorthEast')
    subplot(nf,1,q); q=q+1; hold on;
    plot(TimeAxis, ar,['-' num2str(line(n))]);
    axis tight
    ylabel('Active Rhodopsin');
    subplot(nf,1,q); q=q+1; hold on;
    plot(TimeAxis,p,['-' num2str(line(n))])
    axis tight
    ylabel('pde Activity')
    subplot(nf,1,q); q=q+1; hold on;
    plot(TimeAxis,s,['-' num2str(line(n))])
    axis tight
    ylabel('cGMP synthesis rate')
    subplot(nf,1,q); q=q+1; hold on;
    plot(TimeAxis,g,['-' num2str(line(n))])
    axis tight
    xlabel('Time (sec)')
    ylabel('cGMP')
    subplot(nf,1,q);q=q+1; hold on;
    plot(TimeAxis,c,['-' num2str(line(n))])
    axis tight
    ylabel('[calcium]');
    subplot(nf,1,q);q=q+1; hold on;
    plot(TimeAxis,Cones.Responses(n,:),['-' num2str(line(n))]);
    axis tight
    xlabel('Time (sec)')
    ylabel('Current (pA)')
    
end


%% Part IV: Build Retinal Ganglion Cell

% Here, we wire many photoreceptors together to form a ganglion cell, which
% has both an 'on region' and and 'off region' - that is, a region in which
% light will excite the cell, and promote the firing of action potentials,
% and a region in which light will inhibit the cell, inhibbiting it from
% firing action potentials.  The relative strength of each of these areas
% can be changed either by changing the amount of cones that act as input,
% or adding a coefficient that weighs a given input more of less than other
% inputs.  The example cell below has equal numbers of cells that act as
% inputs to the on-region and off-region, and each has a coefficient of 1.
% Note that the density of cones does not have to be as high as the density
% of cones at the retina, since not every cone in a region will be wired to
% a given ganglion cell.  However, it is well established that neighboring
% cones can have an excitatory or inhibitory effect on each other; this
% complication is not built into my model.

% I also want to point out that I'm using the ambiguous term 'excitatory
% and inhibitory input'.  This is because photoreceptors do not fire action
% potentials, but instead promote the generation of action potentials
% downstream by changing the local extracellular voltage.  There is
% probably some well understood relationship between the output current of
% a photoreceptor and the action potentials of adjacent cells, but I have
% not added this complication into my model.  Instead, I will assume that
% the spike rate of downstream neurons changes linearly with output of a
% photoreceptor.  Since it will be convenient to have a V1 cell with some
% baseline firing rate, I am arbitrily assigning the coefficient of
% proportionality as 1.  Thus, the spike rate of a given ganglion cell will
% be a mean of the input from the photoreceptors.

%%%% User-Defined Parameters %%%%

% Choose diameter of ganglion cell
RFDiameter = .05; % in mm

% Choose number of cones hooked up to a ganglion cell (for now, these will
% be evenly distributed between center and surround, and between L and M cells)
ncones = 64;

% Choose the number of radii on which cones can fall (half of these will be
% center, half surroun - for ease of calculations, stick with even numbers for now)
nradii = 6;


%%%% Calculations %%%%

% Set up individual ganglion cells
bins = linspace(0,RFDiameter/2,nradii*2+1);
radii = bins(2:2:end);
boarders = [bins(nradii+1) bins(end)];
Ganglion.Boarders = boarders;
thetas = linspace(0,2*pi,ncones/2+1);

figure(4); clf;
subplot(1,2,1); hold on;
axis tight
title('Ganglion Receptive Field')
xlabel('Retina (mm)')
ylabel('Retina (mm)')

% Red Cones
[X Y] = meshgrid(thetas(1:2:end),radii);
[centersx centersy] = pol2cart(X(:),Y(:));

for n = 1:numel(centersx)
    [x1 y1] = pol2cart(linspace(0,2*pi,20),sqrt(Cones.CollectingArea)/4);
    x = x1 + centersx(n);
    y = y1 + centersy(n);
    plot(x,y,'r','LineWidth',3)
    [x2 y2] = pol2cart(linspace(0,2*pi,20),Cones.Size/2);
    x = x2 + centersx(n);
    y = y2 + centersy(n);
    plot(x,y,'r')
end

% Green Cones
[X Y] = meshgrid(thetas(2:2:end),radii);
[centersx centersy] = pol2cart(X(:),Y(:));

for n = 1:numel(centersx)
    [x1 y1] = pol2cart(linspace(0,2*pi,20),sqrt(Cones.CollectingArea)/4);
    x = x1 + centersx(n);
    y = y1 + centersy(n);
    plot(x,y,'g','LineWidth',3)
    [x2 y2] = pol2cart(linspace(0,2*pi,20),Cones.Size/2);
    x = x2 + centersx(n);
    y = y2 + centersy(n);
    plot(x,y,'g')
end


% Plot boarders
for n = 1:numel(boarders)
    [bx by] = pol2cart(linspace(0,2*pi,50),boarders(n));
    hold on;
    plot(bx,by,'k--')
    axis equal tight
end



%% Part V: V1 Receptive Field

% Here, I wire several ganglion cells together in a row to create a V1 
% receptive field.  In reality, receptive field inputs are not arranged so 
% neatly, but this is a very simple model.  The stimuli used to test this
% cell will vary only along the x-dimension, and not along the y-dimension.
% Thus, the number of ganglion cells wired together is fairly arbitrary,
% but this model could be extended to build receptive fields that vary in
% both size and location.  Again, the firing rate of this V1 cell will be
% linearly related the firing rates of the input ganglion cells.

%%%% User-Defined Parameters %%%%

% Number of ganglion in a straight line
nganglion = 6;


%%%% Calculations %%%%

% Calculate center of each ganglion (along y axis)
gangcenters = 0:RFDiameter/2:RFDiameter/2*(nganglion-1);
gangcenters = gangcenters - (RFDiameter/2*(nganglion-1))/2;

figure(5); hold on;
%subplot(1,2,2); hold on;
axis tight
title('V1 Receptive Field')
xlabel('Retina (mm)')
ylabel('Retina (mm)')

% L Cones (Center)
Cones.Coordinates.LCones = [];
[X Y] = meshgrid(thetas(3:2:end),radii(1:nradii/2));
[centersx centersy] = pol2cart(X(:),Y(:));

% Duplicate single ganglion into a string of ganglia
centersy = repmat(centersy,numel(gangcenters),1)...
    + reshape(repmat(gangcenters,numel(centersy),1),[],1);
centersx = repmat(centersx,numel(gangcenters),1);
Cones.Coordinates.LCones = cat(1,Cones.Coordinates.LCones,[centersx centersy]);
Ganglion.Center.Coordinates.LCones = [centersx centersy];

% L Cones (Surround)
[X Y] = meshgrid(thetas(3:2:end),radii(nradii/2+1:end));
[centersx centersy] = pol2cart(X(:),Y(:));

% Duplicate single ganglion into a string of ganglia
centersy = repmat(centersy,numel(gangcenters),1)...
    + reshape(repmat(gangcenters,numel(centersy),1),[],1);
centersx = repmat(centersx,numel(gangcenters),1);
Cones.Coordinates.LCones = cat(1,Cones.Coordinates.LCones,[centersx centersy]);
Ganglion.Surround.Coordinates.LCones = [centersx centersy];

% Plot L Cones
for n = 1:numel(Cones.Coordinates.LCones(:,1))
    [x1 y1] = pol2cart(linspace(0,2*pi,20),sqrt(Cones.CollectingArea)/4);
    x = x1 + Cones.Coordinates.LCones(n,1);
    y = y1 + Cones.Coordinates.LCones(n,2);
    plot(x,y,'r','LineWidth',3)
    [x2 y2] = pol2cart(linspace(0,2*pi,20),Cones.Size/2);
    x = x2 + Cones.Coordinates.LCones(n,1);
    y = y2 + Cones.Coordinates.LCones(n,2);
    plot(x,y,'r')
end


% M Cones (Center)
Cones.Coordinates.MCones = [];
[X Y] = meshgrid(thetas(2:2:end),radii(1:nradii/2));
[centersx centersy] = pol2cart(X(:),Y(:));

% Replicate single ganglion into many
centersy = repmat(centersy,numel(gangcenters),1)...
    + reshape(repmat(gangcenters,numel(centersy),1),[],1);
centersx = repmat(centersx,numel(gangcenters),1);
Cones.Coordinates.MCones = cat(1,Cones.Coordinates.MCones,[centersx centersy]);
Ganglion.Center.Coordinates.MCones = [centersx centersy];

% M Cones (Surround)
[X Y] = meshgrid(thetas(2:2:end),radii(nradii/2+1:end));
[centersx centersy] = pol2cart(X(:),Y(:));

% Replicate single ganglion into many
centersy = repmat(centersy,numel(gangcenters),1)...
    + reshape(repmat(gangcenters,numel(centersy),1),[],1);
centersx = repmat(centersx,numel(gangcenters),1);
Cones.Coordinates.MCones = cat(1,Cones.Coordinates.MCones,[centersx centersy]);
Ganglion.Surround.Coordinates.MCones = [centersx centersy];

% Plot M Cones
for n = 1:numel(Cones.Coordinates.MCones(:,1))
    [x1 y1] = pol2cart(linspace(0,2*pi,20),sqrt(Cones.CollectingArea)/4);
    x = x1 + Cones.Coordinates.MCones(n,1);
    y = y1 + Cones.Coordinates.MCones(n,2);
    plot(x,y,'g','LineWidth',3)
    [x2 y2] = pol2cart(linspace(0,2*pi,20),Cones.Size/2);
    x = x2 + Cones.Coordinates.MCones(n,1);
    y = y2 + Cones.Coordinates.MCones(n,2);
    plot(x,y,'g')
end

% Plot boarders
for q = 1:numel(gangcenters)
    for n = 1:numel(boarders)
        [bx by] = pol2cart(linspace(0,2*pi,50),boarders(n));
        by = by+gangcenters(q);
        hold on;
        plot(bx,by,'k--')
        axis equal
    end
end


%% Part VI: Stationary Gratings (Spatial Frequency Tuning Curves)

% We have now created a model of a simple V1 cell, which should behave 
% roughly linearly.  I now want to see if I can recreate the receptive
% field of this cell by taking the inverse Fourier transform of its spatial
% frequency turning curves, just as Movshon et al did.  This section
% flashes a stationary square grating over the receptive field, varying the
% phase to find the optimal response, and varying the spatial frequency of
% these gratings to build up a tuning curve.

% In the Movshon paper, spatial frequency tuning curves are built up by
% stimulating the cell with square gratings that switch between polarities
% (dark to light, light to dark).  Since the V1 cells that they tested
% often had a low enough baseline firing rate that their output firing rate 
% was half-rectified, they used this tactic to uncover the amount of 
% excitatory input from on-regions (dark to light) and off-regions (light
% to dark).  Instead of reversing polarities, I made the assumption that
% the excitatory input of an off-region to dark bar is equivalent to the 
% excitatory input of an on-region to a light bar.  This is very likely a
% poor model for such responses, but it works well enough for a simple
% model.

%%%% User-Defined Parameters %%%%

% Spatial frequencies to be tested
sf = linspace(1/.0025,1/.1,30); % in cycles/mm

% Number of phase shifts to be tested (to find the optimal stimulus)
nphase = 8;


%%%% Calculations %%%%

% Characterize Responses (excitatory or inhibitory)
Lexcite = (-1*(Cones.Responses(1,:) - Cones.Responses(1,1)))+Cones.Responses(1,1);
Mexcite = (-1*(Cones.Responses(2,:) - Cones.Responses(2,1)))+Cones.Responses(2,1);
Sexcite = (-1*(Cones.Responses(3,:) - Cones.Responses(3,1)))+Cones.Responses(3,1);
Linhib = Cones.Responses(1,:);
Minhib = Cones.Responses(2,:);
Sinhib = Cones.Responses(3,:);

% Choose optimal location for stimulus (to save time)
for n = 1:numel(sf)
    
    % Calculate width bar (1/2 of stimulus period)
    barwidth = 1/sf(n)*.5;% Width of light/dark bars
    
    % Test different phase shifts to find the optimal phase
    for q = 1:nphase
        
        barboarders = (-RFDiameter + q*2*barwidth/nphase):barwidth:RFDiameter;
        if barboarders(1) ~= -RFDiameter
            barboarders = [-RFDiameter barboarders];
        end
        if barboarders(end) ~= RFDiameter
            barboarders = [barboarders RFDiameter];
        end
        tempboarders = barboarders - min(barboarders);
        
        % Center LCones in each bin
        tempLC = Ganglion.Center.Coordinates.LCones(:,1) - min(barboarders);
        tempbins = histc(tempLC,tempboarders);
        LCbins = [tempbins(1:end-2)' tempbins(end-1)+tempbins(end)];
        
        % Surround LCones in each bin
        tempLS = Ganglion.Surround.Coordinates.LCones(:,1) - min(barboarders);
        tempbins = histc(tempLS,tempboarders);
        LSbins = [tempbins(1:end-2)' tempbins(end-1)+tempbins(end)];
        
        % Center MCones in each bin
        tempMC = Ganglion.Center.Coordinates.MCones(:,1) - min(barboarders);
        tempbins = histc(tempMC,tempboarders);
        MCbins = [tempbins(1:end-2)' tempbins(end-1)+tempbins(end)];
        
        % Surround MCones in each bin
        tempMS = Ganglion.Surround.Coordinates.MCones(:,1) - min(barboarders);
        tempbins = histc(tempMS,tempboarders);
        MSbins = [tempbins(1:end-2)' tempbins(end-1)+tempbins(end)];
        
        % Sum excited cones and subtract inhibited cones
        lightbins = zeros(1,numel(tempbins)-1);
        darkbins = zeros(1,numel(tempbins)-1);
        lightbins(1:2:end) = 1;
        darkbins(2:2:end) = 1;
        lightbins = logical(lightbins);
        darkbins = logical(darkbins);
        
        % L Cones
        LCexcite = sum(LCbins(lightbins));
        LSexcite = sum(LSbins(darkbins));
        LCsuppress = sum(LCbins(darkbins));
        LSsuppress = sum(LSbins(lightbins));
        
        % M Cones
        MCexcite = sum(MCbins(lightbins));
        MSexcite = sum(MSbins(darkbins));
        MCsuppress = sum(MCbins(darkbins));
        MSsuppress = sum(MSbins(lightbins));
        
        % Total excitation - total suppression
        temptotalLexcite(q) = LCexcite + LSexcite - LCsuppress - LSsuppress;
        temptotalMexcite(q) = MCexcite + MSexcite - MCsuppress - MSsuppress;
        
    end
    
    % Total excitation - total suppression
    totalLexcite(n) = max(temptotalLexcite);
    totalMexcite(n) = max(temptotalMexcite);
    
    clear temptotalLexcite temptotalMexcite
    
end

totalactivity = totalLexcite+totalMexcite;

figure(6); clf; hold on
%plot(sf,totalLexcite,'r')
%plot(sf,totalMexcite,'g--')
plot(sf,totalactivity,'k*-')
title('Spatial Frequency Tuning Curve')
xlabel('Spatial Frequency (cycles/mm)')
ylabel('activity')


%% Part VII: Fourier Analysis

% Fourier transformation
F = fftshift(ifft(totalactivity));

% Reconstruct spatial axis
dx = 1/(max(sf));
SpaceAxis = ((0:(numel(real(F))-1)) - numel(real(F))/2) .* dx;

% Plot
figure(7); clf; hold on;
title('Fourier Transformation of SF Responses')
plot(SpaceAxis,real(F)./max(real(F)),'b-')
plot(SpaceAxis,imag(F)./max(real(F)),'c--')
plot(SpaceAxis,zeros(1,numel(SpaceAxis)),'k--')
legend('Real Component','Imaginary Component','Location','NorthEast')
xlabel('Retina (mm)')
ylabel('Normalized Response Amplitude')

