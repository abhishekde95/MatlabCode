% A MATLAB script for understanding the link between edge processing and
% spatial integration
% Author - Abhishek De, 2/10
close all; clearvars;
plot_counter = 1;

% Reading a grayscale image 
% I = imread('pout.tif');
images = {'pout.tif'; 'coins.png'};

I_orig = imread('pout.tif');
% I_orig = imread('coins.png');
I_color = imread('peppers.png');
I_orig = rgb2gray(I_color);
I = im2double(I_orig)-0.5;

% Designing a filter
F = zeros(4,2); S1 = F; S2 = F;
S1(:,1) = 1; S2(:,2) = -1;
F = S1+S2;

% Linear filtered image 
I_linearv = conv2(I,F,'same');
I_linearv(I_linearv<0) = 0;
I_linear = I_linearv; % Adding the outputs of the horizontal and vertical filters
% I_linear(I_linear>prctile(I_linear(:),90))= 1;
I_linear = I_linear + 0.5-median(I_linear(:));

% Performing the squaring non-linear filtering operation
I_nonlinS1v = conv2(I,S1,'same'); I_S1 = I_nonlinS1v; I_nonlinS1v(I_nonlinS1v<0) = 0; 
I_nonlinS2v = conv2(I,S2,'same'); I_S2 = I_nonlinS2v; I_nonlinS2v(I_nonlinS2v<0) = 0;
I_nonlin = I_nonlinS1v.^2 + I_nonlinS2v.^2; I_nonlin(I_nonlin<0) = 0;
I_nonlin = I_nonlin./max(I_nonlin(:));
I_nonlin(I_nonlin<0) = 0; % Analogous to rectifying non-linearity
% I_nonlin(I_nonlin>prctile(I_nonlin(:),95))= 1;
I_nonlin = I_nonlin + 0.5-median(I_nonlin(:));

% Performing the square-root non-linear filtering operation
I_nonlin2 = sqrt(I_nonlinS1v) + sqrt(I_nonlinS2v); I_nonlin2(I_nonlin2<0) = 0;
I_nonlin2 = I_nonlin2./max(I_nonlin2(:)); % Analogous to rectifying non-linearity 
I_nonlin2(I_nonlin2<0) = 0;
% I_nonlin2(I_nonlin2>prctile(I_nonlin2(:),95))= 1; 
I_nonlin2 = I_nonlin2 + 0.5-median(I_nonlin2(:));

% I want to separate the grayscale image into low- and high-contrast
Drive_image = (I_nonlinS1v+I_nonlinS2v)/2;
I_locontrast = Drive_image; I_hicontrast = Drive_image;
I_locontrast(I_locontrast>median(I_locontrast(:))) = 0;
I_hicontrast(I_hicontrast<median(I_hicontrast(:))) = 0;

% Plotting the figures 
figure(plot_counter); set(gcf,'Name','Analyses of images: grayscale')
subplot(321); imshow(I_orig); axis square; set(gca,'XTick',[],'YTick',[]); title('original pic');
subplot(322); imshow(im2uint8(I_linear)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title('Linear');
subplot(323); imshow(im2uint8(I_nonlin)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title('Squaring non-linearity');
subplot(324); imshow(im2uint8(I_nonlin2)); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); title('Square-root');
subplot(325); imshow(I_hicontrast); colormap('jet'); axis square; set(gca,'XTick',[],'YTick',[]); title('Hi contrast image');
subplot(326); imshow(I_locontrast); colormap('jet'); axis square; set(gca,'XTick',[],'YTick',[]); title('Lo contrast image');
plot_counter = plot_counter + 1;

%%

% Implementing linear filtering in the RGB domain
I_color_filtered = imfilter(im2double(I_color),F,'conv');
I_color_filtered(I_color_filtered>prctile(I_color_filtered(:),95))= 1;

% Implementing squaring non-linear filtering 
I_color_filteredS1 = imfilter(im2double(I_color),S1,'conv'); I_color_filteredS1(I_color_filteredS1<0) = 0;
I_color_filteredS2 = imfilter(im2double(I_color),S2,'conv'); I_color_filteredS2(I_color_filteredS2<0) = 0;
I_color_filteredNL = I_color_filteredS1.^2 + I_color_filteredS2.^2;
I_color_filteredNL(I_color_filteredNL<0) = 0;
I_color_filteredNL = I_color_filteredNL./max(I_color_filteredNL(:));
I_color_filteredNL(I_color_filteredNL>prctile(I_color_filteredNL(:),99)) = 1;


% Implementing square-root non-linear filtering 
I_color_filteredNL2 = sqrt(I_color_filteredS1) + sqrt(I_color_filteredS2);
I_color_filteredNL2(I_color_filteredNL2<0) = 0;
I_color_filteredNL2 = I_color_filteredNL2./max(I_color_filteredNL2(:));
I_color_filteredNL2(I_color_filteredNL2>prctile(I_color_filteredNL2(:),99)) = 1;
figure(plot_counter); set(gcf,'Name','Analyses of images: RGB');
subplot(221); imshow(I_color); axis square; set(gca,'XTick',[],'YTick',[]); title('original pic');
subplot(222); image(im2uint8(I_color_filtered)); axis square; set(gca,'XTick',[],'YTick',[]); title('Linear');
subplot(223); image(im2uint8(I_color_filteredNL)); axis square; set(gca,'XTick',[],'YTick',[]); title('Squaring non-linearity');
subplot(224); image(im2uint8(I_color_filteredNL2)); axis square; set(gca,'XTick',[],'YTick',[]); title('Square-root');
plot_counter = plot_counter + 1;


% Simulating the responses of the 3 filters 
x = linspace(-1,1,51);
[X,Y] = meshgrid(x,x);
Resp_lin = zeros(size(X));
Resp_nonlin = zeros(size(X));
Resp_nonlin2 = zeros(size(X));
for ii = 1:size(X,1)
    for jj = 1:size(X,2)
        img = zeros(4,2);
        
        % Linear filter
        img(:,1) = X(ii,jj); img(:,2) = Y(ii,jj);
        Resp_lin(ii,jj) = img(:)'*F(:);
        
        % Squaring filter
        tmpS1 = max([0 img(:)'*S1(:)]); tmpS2 = max([0 img(:)'*S2(:)]);
        Resp_nonlin(ii,jj) = tmpS1^2 + tmpS2^2;
        
        % Square-root filter
        tmpS1 = max([0 img(:)'*S1(:)]); tmpS2 = max([0 img(:)'*S2(:)]);
        Resp_nonlin2(ii,jj) = sqrt(tmpS1) + sqrt(tmpS2);
    end
end

figure(plot_counter); set(gcf,'Name','Filter Isoresponse contours'); 
subplot(221); contour(x,x,Resp_lin,'Linewidth',2); xlabel('S1'); ylabel('S2'); axis square; set(gca,'Tickdir','out'); axis ij; title('Linear'); colormap('gray');
subplot(222); contour(x,x,Resp_nonlin,'Linewidth',2); xlabel('S1'); ylabel('S2'); axis square; set(gca,'Tickdir','out'); axis ij; title('Squaring'); colormap('gray');
subplot(223); contour(x,x,Resp_nonlin2,'Linewidth',2); xlabel('S1'); ylabel('S2'); axis square; set(gca,'Tickdir','out'); axis ij; title('Square-root'); colormap('gray');
plot_counter = plot_counter + 1;


%% Implementing a different method of processing colore images
I_color_filtered = imfilter(im2double(I_color),F,'conv');
I_color_filtered = squeeze(I_color_filtered(:,:,1)-I_color_filtered(:,:,2));
I_color_filtered(I_color_filtered<0) = 0;
I_color_filtered = I_color_filtered./(1.5*max(I_color_filtered(:)));
I_color_filtered = I_color_filtered + 0.5-median(I_color_filtered(:));

% Implementing squaring non-linear filtering 
I_color_filteredS1 = imfilter(im2double(I_color),S1,'conv'); I_color_filteredS1 = squeeze(I_color_filteredS1(:,:,1)-I_color_filteredS1(:,:,2)); I_color_filteredS1(I_color_filteredS1<0) = 0;
I_color_filteredS2 = imfilter(im2double(I_color),S2,'conv'); I_color_filteredS2 = squeeze(I_color_filteredS2(:,:,1)-I_color_filteredS2(:,:,2)); I_color_filteredS2(I_color_filteredS2<0) = 0;
I_color_filteredNL = I_color_filteredS1.^2 + I_color_filteredS2.^2;
I_color_filteredNL(I_color_filteredNL<0) = 0;
I_color_filteredNL = I_color_filteredNL./(1.5*max(I_color_filteredNL(:)));
I_color_filteredNL = I_color_filteredNL + 0.5-median(I_color_filteredNL(:));

% Visualizing the filter
A = 0.5*ones(100,100,3);
A(20:80,20:50,1) = 1.0; A(20:80,20:50,2) = 0.0; % Subunit 1-> R-G
A(20:80,51:80,1) = 0.0; A(20:80,51:80,2) = 1.0; % Subunit 1-> G-R

% Plotting the images
figure(plot_counter); set(gcf,'Name','Analyses of images: RGB');
subplot(221); image(A); axis square; set(gca,'XTick',[],'YTick',[]); title('Filter');
subplot(222); imshow(I_color); axis square; set(gca,'XTick',[],'YTick',[]); title('original pic');
subplot(223); imshow(im2uint8(I_color_filtered)); axis square; set(gca,'XTick',[],'YTick',[]); title('Linear'); colormap('gray'); hold off;
subplot(224); imshow(im2uint8(I_color_filteredNL)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; title('Squaring non-linearity');

plot_counter = plot_counter + 1;


%% Trying a different approach
% Data - 07/20
% 1) Convert RGB into LMS cone-excitations 2) Normalize cone-signals seprately 
% 3) Then run the spatial filter 

close all; clearvars;
I_color = imread('peppers.png');

% Designing a filter
F = zeros(4,2); S1 = F; S2 = F;
S1(:,1) = 1; S2(:,2) = -1;
F = S1+S2;

% Forming the M matrix 
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

row  = size(I_color,1); col = size(I_color,2);
I_cone_excitation = M*RGB2XWFormat(im2double(I_color))';
I_cone_excitation = XW2RGBFormat(I_cone_excitation',row,col);
I_Lcone_excitation = I_cone_excitation(:,:,1);
I_Mcone_excitation = I_cone_excitation(:,:,2);
I_Scone_excitation = I_cone_excitation(:,:,3);

I_Lcone = (I_Lcone_excitation-median(I_Lcone_excitation(:)))/(median(I_Lcone_excitation(:)));
I_Mcone = (I_Mcone_excitation-median(I_Mcone_excitation(:)))/(median(I_Mcone_excitation(:)));
I_Scone = (I_Scone_excitation-median(I_Scone_excitation(:)))/(median(I_Scone_excitation(:)));
I_LminusM = I_Lcone - I_Mcone;  
I_MminusL = I_Mcone - I_Lcone; 

% Linear filtering
I_color_filtered = imfilter(I_LminusM,F,'conv');
I_color_filtered(I_color_filtered<0) = 0;
I_color_filtered = I_color_filtered./(1.5*max(I_color_filtered(:)));
I_color_filtered = I_color_filtered + 0.5-median(I_color_filtered(:));

% Nonlinear filtering 
I_LminusM(I_LminusM<0) = 0; I_MminusL(I_MminusL<0) = 0;
I_color_filteredNL = I_LminusM.^2 + I_MminusL.^2; 
I_color_filteredNL(I_color_filteredNL<0) = 0;
I_color_filteredNL = I_color_filteredNL./(1.5*max(I_color_filteredNL(:)));
I_color_filteredNL = I_color_filteredNL + 0.5-median(I_color_filteredNL(:));

% Visualizing the filter
A = 0.5*ones(100,100,3);
RGB_forLminusM = inv(M')*[1;-1;0]; RGB_forLminusM = RGB_forLminusM./(max(abs(RGB_forLminusM)));
RGB_forMminusL = inv(M')*[-1;1;0]; RGB_forMminusL = RGB_forMminusL./(max(abs(RGB_forMminusL)));
A(20:80,20:50,1) = RGB_forMminusL(1); A(20:80,20:50,2) = RGB_forMminusL(2); A(20:80,20:50,3) = RGB_forMminusL(3); % Subunit 1-> L-M
A(20:80,51:80,1) = RGB_forLminusM(1); A(20:80,51:80,2) = RGB_forLminusM(2); A(20:80,51:80,3) = RGB_forLminusM(3); % Subunit 2-> M-L

% Normalizing the L-M  and M-L responses
I_LminusM_norm = I_LminusM - median(I_LminusM(:));
I_LminusM_norm = I_LminusM_norm./(2*max(I_LminusM_norm(:)));
I_LminusM_norm = I_LminusM_norm + 0.5;

I_MminusL_norm = I_MminusL - median(I_MminusL(:));
I_MminusL_norm = I_MminusL./(2*max(I_MminusL(:)));
I_MminusL_norm = I_MminusL_norm + 0.5;


% Plotting the images
plot_counter = 1;
figure(plot_counter); set(gcf,'Name','Analyses of images: RGB');
subplot(321); imagesc(A); axis square; set(gca,'XTick',[],'YTick',[]); title('Filter');
subplot(322); imshow(I_color); axis square; set(gca,'XTick',[],'YTick',[]); title('original pic');
subplot(323); imshow(im2double(I_LminusM_norm)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); title('Rectified L-M map');
subplot(324); imshow(im2double(I_MminusL_norm)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); title('Rectified M-L map');
subplot(325); imshow(im2double(I_color_filtered)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; title('Linear filter');
subplot(326); imshow(im2double(I_color_filteredNL)); axis square; set(gca,'XTick',[],'YTick',[]); colormap('gray'); hold off; title('Squaring non-linearity');
plot_counter = plot_counter + 1;

% Trying out the non-rectified maps
I_LminusMnon = I_Lcone - I_Mcone;  
I_MminusLnon = I_Mcone - I_Lcone;
I_LmminusMnon = I_LminusMnon - median(I_LminusMnon(:));
I_LminusMnon = I_LminusMnon./(2*max(abs(I_LminusMnon(:))));
I_LminusMnon = I_LminusMnon + 0.5;

I_MmminusLnon = I_MminusLnon - median(I_MminusLnon(:));
I_MminusLnon = I_MminusLnon./(2*max(abs(I_MminusLnon(:))));
I_MminusLnon = I_MminusLnon + 0.5;

% PLotting the distribution of the L-, M- and S-cone contrast
figure(plot_counter);
subplot(221); histogram(I_Lcone_excitation(:),0:0.01:0.2,'DisplayStyle','stairs','EdgeColor',[1 0 0]); hold on;
histogram(I_Mcone_excitation(:),0:0.01:0.2,'DisplayStyle','stairs','EdgeColor',[0 1 0]); 
histogram(I_Scone_excitation(:),0:0.01:0.2,'DisplayStyle','stairs','EdgeColor',[0 0 1]); 
plot(median(I_Lcone_excitation(:)),0,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Mcone_excitation(:)),0,'v','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Scone_excitation(:)),0,'v','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('cone excitation'); ylabel('Pixel count'); legend('L','M','S');
subplot(222); histogram(I_Lcone(:),-1:0.1:1,'DisplayStyle','stairs','EdgeColor',[1 0 0]); hold on;
histogram(I_Mcone(:),-1:0.1:1,'DisplayStyle','stairs','EdgeColor',[0 1 0]); 
histogram(I_Scone(:),-1:0.1:1,'DisplayStyle','stairs','EdgeColor',[0 0 1]); 
plot(median(I_Lcone(:)),0,'v','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Mcone(:)),2500,'v','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
plot(median(I_Scone(:)),5000,'v','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out'); xlabel('cone contrast'); ylabel('Pixel count'); legend('L','M','S');
subplot(223); imshow(im2double(I_LminusMnon)); axis square; colormap('gray'); set(gca,'XTick',[],'YTick',[]); title('L-M map');
subplot(224); imshow(im2double(I_MminusLnon)); axis square; colormap('gray'); set(gca,'XTick',[],'YTick',[]); title('M-L map');

plot_counter = plot_counter + 1;





