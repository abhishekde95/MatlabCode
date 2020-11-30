function cs_density=VisAngle2Retina_Ecc(Eccentricity)
% function cs_density=VisAngle2Retina_Ecc(Eccentricity)
% Calculates cone density given eccentricity for macaque retina in
% cones/mm^2
% Taken from Goodchild et al., 1996
coeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; %for macaque (Goodchild et al., 1996)
cs_density= (coeffs(1)*(exp(coeffs(2)*Eccentricity)))+...
            (coeffs(3)*(exp(coeffs(4)*Eccentricity)))+...
            (coeffs(5)*(exp(coeffs(6)*Eccentricity)));
cs_density=cs_density*1e3; %in cones/mm^2 Axis Label left out in Goodchild et al., 1996?
end