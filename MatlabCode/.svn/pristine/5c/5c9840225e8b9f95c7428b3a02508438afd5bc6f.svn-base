function [return_basis_vec, plot_counter, subunits,sub_rgb_idx] = compute_basis_vec(plot_counter,nrandnums_perchannel,mask)
% Computes and plots the STAs, PCs and the temporal variation of the
% components of the basis vector.

global nstixperside maxT

out = STCOVmex('return'); % returns the covariance matrix on frame by frame basis
STS = out{1};  % A (dimension) x 9(frames) matrix
STCross = out{2};  % A (dimension).^2 x 9(frames) matrix
nspikes = out{3}; % Number of spikes in the given file
clear STCOVmex;
clear out;

% Coverting the STS and the STCross into STA and STC respectively
STAs = STS/nspikes;
STCs = zeros(size(STCross));
for i = 1:maxT
    tmp = STS(:,i)*STS(:,i)';
    STCs(:,i) = (nspikes.*STCross(:,i)-tmp(:))/(nspikes*(nspikes-1));
end

an_mask = [];
an_mask = zeros(nrandnums_perchannel, 1);   % Pixel mask. Someday make a tool to make this non-zero.
Lmask = logical(repmat(~an_mask(:),[3 1]));
PCs = [];
eig_PC = zeros(3,size(STCs,2));   % to store the eigenvalues of  the first 3 principle components across all frames

% Obtaining the eigenvectors and their corresponding eigenvalues
for i = 1:size(STCs,2)
    STC = reshape(STCs(:,i), 3*nrandnums_perchannel, 3*nrandnums_perchannel);
    subSTC = STC(Lmask, Lmask);
    subSTA = STAs(Lmask,i);
    
    % subtracting the samples from the STA to ensure the PCs are
    % orthogonal to the the STA.
    P = eye(size(subSTC)) - subSTA*inv(subSTA'*subSTA)*subSTA';
    subSTC = P*subSTC*P';
    [tmp,d] = eig(subSTC);
    v = repmat(double(Lmask),[1 size(tmp,2)]);
    v(Lmask,:) = real(tmp);
    [~, idxs] = sort(diag(d));
    v = v(:,idxs);
    v = v(:,[end end-1 end-2]);  % Collecting the first 3 principle components
    PCs = cat(3, PCs, v);
    f = sort(diag(d));
    eig_PC(:,i) = [f(end);f(end-1);f(end-2)];
end
save('eig_PC.mat','eig_PC');

% flipping the signs of the principal components to make them consistent
for i = size(PCs,2):-1:1
    for j = 1:size(PCs,3)-1
        PCs(:,i,j+1) = sign(squeeze(PCs(:,i,j))'*squeeze(PCs(:,i,j+1)))* PCs(:,i,j+1);
    end
end

% Displaying the temporal profile across the frames of the eigenvalues for
% the first 3 principle components
figure(plot_counter);
plot(eig_PC(1,:)/max(eig_PC(1,:)),'Linewidth',2); hold on;
plot(eig_PC(2,:)/max(eig_PC(2,:)),'r','Linewidth',2); hold on;
plot(eig_PC(3,:)/max(eig_PC(3,:)),'g','Linewidth',2); hold on;
xlabel('frame number'); ylabel('Normalized PC1 eigenvalue');
legend('PC 1','PC 2','PC 3');
title('Temporal variation of PC1 eigenvalue');
plot_counter = plot_counter + 1;

muvect = repmat([.5 ;.5; .5],nstixperside^2,1); % creating a 300 x 1 array with each entry as 0.5
% The plotting begins here

FigHandle = figure(plot_counter);
set(FigHandle, 'Position', [50, 70, 1449, 705]);
for i = 1:size(STAs,2) % evaluates for each frame
    
    normfactor = 0.5/((max(abs(STAs(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
    STA = STAs(:,i);
    STA = expand_vector(STA,nrandnums_perchannel,mask,1);
    STA = normfactor*(STA)+muvect;  % This makes the values fall back within a range of 0 and 1.
    STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
    subplot(5,size(STAs,2),i);
    image(STA); % for viewing the image
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (i == 1)
        ylabel('STA');
    end
    for j = 1:size(v,2) % size(v,2) denotes the number of principle components we are interested in looking at
        PC_int = PCs(:,j,i);
        PC_int = 0.5*(PC_int)/(max(abs(PC_int))*1.05);
        PC = expand_vector(PC_int,nrandnums_perchannel,mask,1);
        normfactor = 1;
        PC = normfactor*PC+muvect;
        PC = reshape(PC,[nstixperside nstixperside 3]);
        subplot(5,size(STAs,2), j*size(STAs,2)+i);
        image(PC);
        if (i == 1)
            ylabel(strcat('PC',num2str(j)));
        end
        set(gca,'XTick',[],'YTick',[]); axis square;
    end
    
    STC = reshape(STCs(:,i),[sqrt(length(STCs(:,i))),sqrt(length(STCs(:,i)))]);
    STV = diag(STC);
    STV = 0.5*(STV)/(max(abs(STV))*1.05);
    STV = expand_vector(STV,nrandnums_perchannel,mask,1);
    normfactor = 1;
    STV = normfactor*STV+muvect;
    STV = reshape(STV,[nstixperside nstixperside 3]);
    subplot(6,size(STAs,2),5*size(STAs,2)+i);
    image(STV);
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (i == 1)
        ylabel('STV');
    end
end
plot_counter = plot_counter + 1;

%**************************************************************************
% plotting the temporal variation of RGB components of STA of each subunits
%**************************************************************************
subunits = nrandnums_perchannel;
n_elements = nstixperside^2;
sub_rgb_idx = cell(3,subunits); % for 3 color channels (RGB)
temp_plot = zeros(3,maxT,subunits);

% The return_basis_vec is the variable where u feed in the basis vector of your interest.
% A word of caution. The basis vector should be of size 6 x 9. Below is a list
% of choices. Comment the basis vectors which u are not going to use.
% Defining basis vector here

prompt = ['Choose a basis vector: (1,2,3,4) \n'...
    ' 1) STA\n '...
    '2) Principal Component 1\n '...
    '3) Principal Component 2\n '...
    '4) Principal Component 3\n '];
str = input(prompt,'s');
choice = str2num(str);

if (choice == 1)
    return_basis_vec = STAs;
else
    return_basis_vec = squeeze(PCs(:,choice-1,:));
end


% A different analysis starts here
frame_of_interest = 5;
prompt = 'Do you want to modify the current basis vector? (Y/N)\n';
str = input(prompt,'s');
if ((strcmp('Y',str) || strcmp('y',str)))
    prompt1 = 'Would u like to try the SVD version of separating color and spatial basis vector? (Y/N)\n';
    str1 = input(prompt1,'s');
end

% basis_vec = squeeze(PCs(:,3,:));
figure(plot_counter);
sub_rgb_idx = [];
for i = 1:subunits
    sub_rgb_idx = [sub_rgb_idx; i, subunits+i, 2*subunits+i]; %R,G,B
end

if ((strcmp('Y',str) || strcmp('y',str)))
    % Enter only if u want to modify the return_basis_vec
    if (strcmp('Y',str1) || strcmp('y',str1))
        % using SVD to extract a set of separate chromatic and spatio-temporal basis functions
        % CODE BROKEN - ABHISHEK HACK HERE, Need to define basis vector here
               
        [U,S,V] = svd(return_basis_vec);
        % U is a 6x6 matrix whereas V is a 9x9 matrix. The 1st column
        % of U is the spatio-temporal basis function (highest eigenvalue)
        % and 1st row of V is the contrast basis function (highest
        % eigenvalue)
        prompt2 = 'Choose a SVD principle component. (1,2,3)\n';
        comp_num = str2num(input(prompt2,'s'));        
        contrast_vec = U(:,comp_num); % chromatic receptive field
        space_temp_vec = V(:,comp_num); % spatiotemporal receptive field
        
        for i = 1:subunits
            for j = 1:size(return_basis_vec,2)
                return_basis_vec(sub_rgb_idx(i,1),j) = space_temp_vec(j) * contrast_vec(sub_rgb_idx(i,1)); % R plane
                return_basis_vec(sub_rgb_idx(i,2),j) = space_temp_vec(j) * contrast_vec(sub_rgb_idx(i,2)); % G plane
                return_basis_vec(sub_rgb_idx(i,3),j) = space_temp_vec(j) * contrast_vec(sub_rgb_idx(i,3)); % B plane
                temp_plot(1,j,i) = return_basis_vec(sub_rgb_idx(i,1),j);
                temp_plot(2,j,i) = return_basis_vec(sub_rgb_idx(i,2),j);
                temp_plot(3,j,i) = return_basis_vec(sub_rgb_idx(i,3),j);
            end
            subplot(subunits,1,i); % plotting the R:G, G:B and B:R ratio
            plot(temp_plot(1,:,i),'r'); hold on;
            plot(temp_plot(2,:,i),'g'); hold on;
            plot(temp_plot(3,:,i)); hold on;
            xlabel('frames'),ylabel('Intensity');title('Time Course of Basis Vectors'); hold on;
        end
        hold off;
        fprintf('Index of separability is %d\n\n', S(1,1)/sum(S(:)));
        
    else
        % This is where u introduce a new basis function where the RGB
        % ratio is constant across the time frames
        for i = 1:subunits
            r = return_basis_vec(sub_rgb_idx(1),frame_of_interest);
            g = return_basis_vec(sub_rgb_idx(2),frame_of_interest);
            b = return_basis_vec(sub_rgb_idx(3),frame_of_interest);
            color = [r g b];
            [max_val,max_ind] = max(abs(color));
            if (max_val ~= color(max_ind))
                max_val = max_val * -1;
            end
            r = r/max_val;
            g = g/max_val;
            b = b/max_val;
            for j = 1:size(return_basis_vec,2)
                temp_basis_vec = return_basis_vec(:,j);
                temp_basis_vec(sub_rgb_idx(1)) = temp_basis_vec(sub_rgb_idx(1)) * r;
                temp_basis_vec(sub_rgb_idx(2)) = temp_basis_vec(sub_rgb_idx(2)) * g;
                temp_basis_vec(sub_rgb_idx(3)) = temp_basis_vec(sub_rgb_idx(3)) * b;
                return_basis_vec(:,j) = temp_basis_vec;
            end
            fprintf('Introducing a new basis vector for a luminance cell for subunit %d \n', i);
            for j = 1:size(return_basis_vec,2)
                temp_plot(1,j,i) = return_basis_vec(sub_rgb_idx(1),j)/return_basis_vec(sub_rgb_idx(2),j); % R:G
                temp_plot(2,j,i) = return_basis_vec(sub_rgb_idx(2),j)/return_basis_vec(sub_rgb_idx(3),j); % G:B
                temp_plot(3,j,i) = return_basis_vec(sub_rgb_idx(3),j)/return_basis_vec(sub_rgb_idx(1),j); % B:R
            end
            subplot(subunits,1,i);
            plot(temp_plot(1,:,i),'r'); hold on;
            plot(temp_plot(2,:,i),'g'); hold on;
            plot(temp_plot(3,:,i)); hold on;
            xlabel('frames'),ylabel('Intensity');legend('R:G','G:B','B:R'); hold on;
        end
    end
else
    for i = 1:subunits
        temp_plot(1,:,i) = return_basis_vec(i,:); % R
        temp_plot(2,:,i) = return_basis_vec(subunits+i,:); % G
        temp_plot(3,:,i) = return_basis_vec(2*subunits+i,:); % B
        
        subplot(subunits,1,i);
        plot(temp_plot(1,:,i),'r','LineWidth',2); hold on;
        plot(temp_plot(2,:,i),'g','LineWidth',2); hold on;
        plot(temp_plot(3,:,i),'LineWidth',2); hold on;
        xlabel('frames'),ylabel('Intensity');legend('R','G','B'); title('Time Course of Basis Vectors'); hold on;
    end
end

hold off;
plot_counter = plot_counter + 1;

% Displaying the basis vectors
FigHandle = figure(plot_counter);
set(FigHandle, 'Position', [50, 250, 1449, 350]);
for i = 1:maxT % evaluates for each frame
    normfactor = 0.5/((max(abs(return_basis_vec(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5
    tmp = return_basis_vec(:,i);
    tmp = expand_vector(tmp,nrandnums_perchannel,mask,1);
    tmp = normfactor*(tmp)+muvect; % This makes the values fall back within a range of 0 and 1.
    tmp = reshape(tmp,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
    subplot(1,maxT,i);
    image(abs(tmp)); % for viewing the image
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (i == 1)
        ylabel('Basis Functions');
    end
end
plot_counter = plot_counter + 1;

end

