function [val] = gauss_fun(x,E,data)

S = x(1); %width
D = x(2);   %shift
A = x(3);   %peak
V = x(4); %vertical shift

val = 0;
for i = 1:length(E)
    val = val + ((((A*exp(-(E(i)-D).^2/(2*S.^2)) / (S*sqrt(2*pi)))+V)) - data(i))^2;
end
