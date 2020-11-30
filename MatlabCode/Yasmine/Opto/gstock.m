function gstock(N)


%% Yasmine El-Shamayleh
% Locate my glycerol stock!

%% WHERE'S MY STOCK

close all;
figure;

dim = 9; 
BOX = zeros(dim,dim);

BOX(N) = 1;
BOX = BOX';

%% PLOT

imagesc(BOX); axis square;
set(gca,'tickDir','out');
set(gca,'xtick',1:9,'ytick',1:9);
xlabel('COLUMN#');
ylabel('ROW#');

text =[];
text{1} = sprintf('STOCK #  %.0f ', N);
% text{2} = sprintf('ROW #  .0%f', R);
% text{3} = sprintf('COLUMN #  .0%f', C);

title(text);

