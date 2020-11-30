function [LMS1, LMS2] = getLMSvals(net_spectra1, net_spectra2, fundamentals)
L1 = net_spectra1 * fundamentals(:,1);
M1 = net_spectra1 * fundamentals(:,2);
S1 = net_spectra1 * fundamentals(:,3);
L2 = net_spectra2 * fundamentals(:,1);
M2 = net_spectra2 * fundamentals(:,2);
S2 = net_spectra2 * fundamentals(:,3);
LMS1 = [L1; M1; S1];
LMS2 = [L2; M2; S2];
end

