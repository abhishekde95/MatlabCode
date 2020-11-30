function [fit] = lorentzsum_poles(beta, x)

fit = abs(beta(1)) ./ (1 + (x ./ abs(beta(2))).^2).^beta(3);
fit = (fit + abs(beta(4)) ./ (1 + (x ./ abs(beta(5))).^beta(6)));
