function [response] = neuron(currentLMS,noise,stimIntensity)
global ConeWeights


if nargin < 3
    stimIntensity = 1;
end


%% Linear Functions

% response = max((ConeWeights(1:2) * currentLMS'),3);

% response = (stimIntensity*(ConeWeights(1:2)*currentLMS'))+10;


%% Continuous Non-Linear Functions

% Elipse
 response = (currentLMS(:,1).^2*ConeWeights(1).^2 + currentLMS(:,2).^2*ConeWeights(2).^2);
 
% response = ConeWeights(1) + ConeWeights(2)*currentLMS(1) + ConeWeights(3)*currentLMS(2) + ConeWeights(4)*currentLMS(1).^2 ...
%    + ConeWeights(5)*currentLMS(2).^2 + ConeWeights(6)*currentLMS(1)*currentLMS(2);

% response = abs((ConeWeights(1:2) * currentLMS').^3);

%  response = (ConeWeights(1:2) * currentLMS') .* 10;
%  response(ConeWeights(1:2) * currentLMS' < 0) = 0;
%  response(ConeWeights(1:2) * currentLMS' > .8) = 8;


% if currentLMS(1) < 0
%     response = 0;
% else
%     response = currentLMS(1);
% end

% response = (max(0,ConeWeights(1:2) * currentLMS')).^2;

% if ConeWeights * currentLMS' < 0
%     response = (ConeWeights * currentLMS').^2;
% else
%     response = abs(ConeWeights * currentLMS').^3;
% end

% response = (ConeWeights(1:2) * currentLMS').^2;

% response = abs(ConeWeights(1:2) * currentLMS').^3;


%% Discontinuous Non-Linear Functions

% response = abs(2*(currentLMS(1)).^2) + currentLMS(2);

% if currentLMS(1) < .5
%     response = 0;
% else
%     response = currentLMS(1).^4 + currentLMS(2).^4 + currentLMS(3).^4;
% end

% if sum(currentLMS > .3 & currentLMS < .7) == 3
%     response = 1;
% else
%     response = 0;
% end

%% Addative noise option

if noise == 1
    response = response + normrnd(0,1,size(response));
elseif noise == 2
    response = poissrnd(max(0,response));
elseif noise == 3
    response = response + normrnd(0,2*abs(response),size(response));
end


