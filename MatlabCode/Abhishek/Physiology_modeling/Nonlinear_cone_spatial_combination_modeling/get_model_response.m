function [Out,driveS1,driveS2,net_drive] = get_model_response(Inp,RF,S1,S2,cone_signal_combination)
max_FR = 100;
% keyboard;
InpS1 = Inp(S1,:); InpS2 = Inp(S2,:);
RFS1 = RF(S1,:); RFS2 = RF(S2,:);
if cone_signal_combination == 1
    % Linear summation of cone signals
    driveS1 = InpS1(:)'*RFS1(:);
    driveS2 = InpS2(:)'*RFS2(:);
    net_drive = driveS1+driveS2;
elseif cone_signal_combination == 2
    % Ellipsoid
    driveS1 = InpS1.*RFS1; driveS2 = InpS2.*RFS2;
    tmp_drive = [driveS1; driveS2];
    net_drive = sum((sum(tmp_drive,1)).^2);
elseif cone_signal_combination == 3
    % Hyperboloid
    driveS1 = InpS1.*RFS1; driveS2 = InpS2.*RFS2;
    tmp_drive = [driveS1; driveS2];
    tmp_drive = (sum(tmp_drive,1)).^2;
    net_drive = tmp_drive(1) - (tmp_drive(2)/3) + tmp_drive(3);
end
FR = max_FR * 1./(1+exp(-1*(net_drive-1))); % Spike generation mechanism: sigmoidal non-linearity 
Out = FR;
end

