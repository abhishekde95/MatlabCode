%% trying to find the monitor gamut in L/M CC space
% Created Oct_2011 by CH and JPW

clear all
close all

load rig1mon
Mmtx = monspd*fundamentals;
bkgndlms = Mmtx*bkgndrgb;


thetas = linspace(0,2*pi, 5000);
gamut = nan(size(thetas));
incSize = 0.005;
for a = 1:numel(thetas);
    guess = incSize;
    while(1)
        LMS = [cos(thetas(a)), sin(thetas(a)), .5] .* guess; %LMS of vector in cone contrast units
        rgb = Mmtx \ (bkgndlms(:).*((LMS(:)+1)));        
        
        if any(rgb<0) || any(rgb>1)
            gamut(a) = guess-incSize;
            break
        else
            guess = guess+incSize;
        end
    
    end
end


figure
plot(thetas*(180/pi), gamut)
xlim([0, 360])

figure
polar(thetas, gamut)
