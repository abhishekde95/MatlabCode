% 5/8/18 not sure when this was written or whether it works, but it seems
% like a stub. Probably not useful, but who knows
%almost iterate and plot files, but now for isosamp packs
%remove first file in the list (that's the isosamp file), save for later
%run iterate and plot files for the single file

%gather and plot isosamp data on same axes as tested data
%isoStro = nex2stro(<insert path of first file>)
isoData = isoStro.trial(:, [20, 21, 17]); %L, M, TF
hold on;
plot3(isoData(:,1), isoData(:,2), isoData(:,3), 'rs');


%to test whether isoData tf was within 1-1.5 of actual threshold tf
%[tempTheta, tempR] = cart2pol(gregData(:,1), gregData(:,2)); %L and M coordinates of the collected data
%tempR = tempR * 1.5;
%[confirmX, confirmY] = pol2cart(tempTheta, tempR);
%hold on;
%plot3(confirmX, confirmY, data(:,3), 'b.') %this is plotting 1.5 times the
%LM coords but keeping the same tf to see if this aligns with what isosamp collected

%Lopp = sign(gregData(:,1))~=sign(gregData(:,2))
[tempTheta, tempR] = cart2pol(gregData(:,1), gregData(:,2));
LoppTheta = tempTheta >= pi/2;
figure
plot(gregData(LoppTheta, 3), tempR(LoppTheta), '.')
hold on;
plot(gregData(~LoppTheta, 3), tempR(~LoppTheta), 's')

%plot LM coords x axes, tf y axes for real data and iso seed data