% This is a pet project to understand if a non linear cell could give rise
% to line weighting functions as described in Movshon et al;1978
% Author - Abhishek De

close all; clearvars;
% The first thing to do would be to create a liner cell/ a gabor filter and
% calculate the lineweighting function
L = 100;
x = linspace(-1,1,L);
[x,y] = meshgrid(x,x);
gamma = 1/1; % Aspect ratio of the gabor;
sigma = 0.5;
theta = pi/2;
phase = pi/2;
sf = 0.5; % in cycles per degree
X = x*cos(-theta) + y*sin(-theta);
Y =-x*sin(-theta) + y*cos(-theta);
expterm = exp(-(X.^2 + gamma^2*Y.^2)/2/sigma^2);
costerm = cos(bsxfun(@plus, 2*pi*Y*sf, phase));
gaborfilter = expterm.*costerm;
% Now marking the two subunits
cutoffvalue = 0.0;
[a,b] = find(gaborfilter>=cutoffvalue);
[c,d] = find(gaborfilter<=-cutoffvalue);
masksub1 = zeros(size(gaborfilter));
masksub2 = masksub1;
masksub1(a,b) = 1;
masksub2(c,d) = 1;
sub1 = masksub1.*gaborfilter;
sub2 = masksub2.*gaborfilter;

% Now running a loop to generate a lineweighting function
output_lin = [];
output_nonlin = [];
barwidth = 5;
pow = 1;
for ii = 1:L/barwidth
    input = zeros(size(gaborfilter));
    input(:,(ii-1)*barwidth+1:ii*barwidth) = 1;
    % Spatial summation by a linear cell
    output_lin = [output_lin; exp(0.01*(input(:)'*gaborfilter(:)))];
    % Spatial summation by a nonlinear cell
    tmp_output = (sum(sum(input.*sub1)))^pow + (sum(sum(input.*sub2)))^pow;
    output_nonlin = [output_nonlin; exp(0.01*tmp_output)];
end
figure(1), subplot(221);imagesc(gaborfilter), set(gca,'XTick',[],'YTick',[]); title('RF')
subplot(222); plot(output_lin,'-o','Linewidth',2); set(gca,'Yscale','log');title('Line weighting function');
subplot(223); imagesc(sub1); set(gca,'XTick',[],'YTick',[]);
subplot(224); plot(output_nonlin,'-go','Linewidth',2); set(gca,'Yscale','log');