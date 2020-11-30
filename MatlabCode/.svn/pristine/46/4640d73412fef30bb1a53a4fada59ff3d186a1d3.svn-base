% Short script for getting the STAs for color cells
% Author - Abhishek De, 1/17

close all; clearvars;
load Output_List_Gun.mat
R = size(Output_List,1);
num_rows = 10; % Number of cells in a figure
nstixperside = 10;
idxs = [1;2;3;4;6;7;10;12;13;14;15;20;21;24;26;29;30;31;33;34;37;38;40;41;...
    46;48;50;51;53;56;57;59;61;63;65;72;73;76;77;78;81;84;88;90;92;95;96;99;100];
idxs2 = [17;24;39;6;7;8;43;36;45];
figure(1);
for jj = 1:numel(idxs2)
    % gunnoise based STA
    ii = idxs(idxs2(jj));
    tmp_vec_gun = Output_List{ii,3};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*Output_List{ii,3} + 0.5;
    im = reshape(im,[nstixperside nstixperside 3]);
    subplot(3,3,jj); image(im); set(gca,'XTick',[],'YTick',[]);
end