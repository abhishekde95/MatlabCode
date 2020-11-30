function [Out,driveS1,driveS2,net_drive] = get_model_spatial_response(Inp,RF,S1,S2,spatial_summation)
max_FR = 100;
InpS1 = Inp(S1,:); InpS2 = Inp(S2,:);
RFS1 = RF(S1,:); RFS2 = RF(S2,:);
if spatial_summation == 1
    % Linear 
    driveS1 = InpS1(:)'*RFS1(:);
    driveS2 = InpS2(:)'*RFS2(:);
    net_drive = driveS1+driveS2;
elseif spatial_summation == 2
    % Squaring and adding 
    driveS1 = InpS1(:)'*RFS1(:); driveS2 = InpS2(:)'*RFS2(:);
    net_drive = (driveS1).^2 + (driveS2).^2;
elseif spatial_summation == 3
    % Hyperbolic
    driveS1 = InpS1(:)'*RFS1(:); driveS2 = InpS2(:)'*RFS2(:);
    net_drive = driveS1.*driveS2;
elseif spatial_summation == 4
    % Squaring and subtracting
    driveS1 = InpS1(:)'*RFS1(:); driveS2 = InpS2(:)'*RFS2(:);
    net_drive = (driveS1).^2 - ((driveS2).^2/10);
end
FR = max_FR * 1./(1+exp(-1*(net_drive-1))); % Spike generation mechanism: sigmoidal non-linearity 
Out = FR;
end

