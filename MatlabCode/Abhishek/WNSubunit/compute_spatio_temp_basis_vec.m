function [ch,new_vec,return_basis_vec, plot_counter,subunits,sub_rgb_idx] = compute_spatio_temp_basis_vec(plot_counter,nrandnums_perchannel,mask,flag,mode)

global nstixperside maxT M
% Here dimension = nframes x ncolor x nsubunits

out = STCOV_st('return'); % returns the covariance matrix on frame by frame basis
STS = out{1};  % A (dimension) x 9(frames) matrix
STCross = out{2};  % A (dimension x frames)x (dimension x frames)  matrix
nspikes = out{3}; % Number of spikes in the given file
clear STCOV_st;
clear out;
% Coverting the STS and the STCross into STA and STC respectively
STAs = STS/nspikes;
tmp = STS(:)*STS(:)';
STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));
% keyboard;

% Obtaining the eigenvectors and their corresponding eigenvalues  
% subtracting the samples from the STA to ensure the PCs are
% orthogonal to the the STA.
P = eye(size(STCs)) - STAs(:)*inv(STAs(:)'*STAs(:))*STAs(:)'; % WHAT DOES THIS LINE MEAN
STCs = P*STCs*P';
[tmp,d] = eig(STCs);
eig_PC = sort(diag(d)); % storing all the eigenvalues
figure(plot_counter); % plotting the eigenvalues
plot(flipud(eig_PC),'*'), ylabel('Eigenvalues'),xlabel('index'),title('Principle components');
plot_counter = plot_counter + 1;

v = real(tmp);
[~, idxs] = sort(diag(d));
v = v(:,idxs);
suppresive_PC = 2;
PCs = v(:,[end end-1 end-2 suppresive_PC]);  % Collecting the first 3 principle components and a suppresive PC(Acc to Rust et al; 2005)

% -----------------------------------------------------------------------------------------------------
% Flipping the STAs and the PCs such that the last frame appears first and the first frame appears last
tmp1 = [];
for i=1:size(PCs,2)
    tmp2 = fliplr(reshape(PCs(:,i),[nrandnums_perchannel*3 maxT]));
    tmp1 = [tmp1, tmp2(:)];
end
PCs = tmp1;
STAs = fliplr(STAs); 
% -----------------------------------------------------------------------------------------------------

muvect = repmat([.5 ;.5; .5],nstixperside^2,1); % creating a 300 x 1 array with each entry as 0.5
% The plotting begins here


FigHandle = figure(plot_counter);
set(FigHandle, 'Position', [50, 70, 1449, 705]);
for i = 1:size(STAs,2) % evaluates for each frame
    normfactor = 0.5/((max(abs(STAs(:)))) + 0.0001); % This step is necessary to constrain the values within [-0.5, 0.5]
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
end
for i = 1:size(PCs,2) % size(v,2) denotes the number of principle components we are interested in looking at
    PC_int = PCs(:,i);
    PC_int = 0.5*(PC_int)/(max(abs(PC_int))*1.05);
    PC_int = reshape(PC_int,[nrandnums_perchannel*3 maxT]);
    for j = 1:maxT
        PC = expand_vector(PC_int(:,j),nrandnums_perchannel,mask,1);
        normfactor = 1;
        PC = normfactor*PC+muvect;
        PC = reshape(PC,[nstixperside nstixperside 3]);
        subplot(5,maxT, i*maxT+j);
        image(PC);
        if (j == 1 && i < size(PCs,2))
            ylabel(strcat('PC',num2str(i)));
        elseif (j == 1 && i == size(PCs,2))
            ylabel(strcat('PC:end-',num2str(suppresive_PC-1)));
        end
        set(gca,'XTick',[],'YTick',[]); axis square;
    end
end
plot_counter = plot_counter + 1;

%**************************************************************************
% plotting the temporal variation of RGB components of STA of each subunits
%**************************************************************************
subunits = nrandnums_perchannel;
n_elements = nstixperside^2;
sub_rgb_idx = cell(3,subunits); % declaring a cell for 3 color channels (RGB)
temp_plot = zeros(3,maxT,subunits);
LMS_temp_plot = zeros(3,maxT,subunits);

% The return_basis_vec is the variable where u feed in the basis vector of your interest.
% A word of caution. The basis vector should be of size 6 x 9. Below is a list
% of choices. Comment the basis vectors which u are not going to use.
% Defining basis vector here
new_vec = []; ch = [];
if flag == 1
    return_basis_vec = STAs; % for use onlyu within this function
    prompt = 'Select the first vector. (STA->1, PC1->2, PC2->3, PC3->4, PC_supp->5)';
    ch1 = input(prompt,'s');
    prompt = 'Select the second vector. (STA->1, PC1->2, PC2->3, PC3->3, PC_supp->5)';
    ch2 = input(prompt,'s');
    ch = [str2num(ch1), str2num(ch2)];
    tmp1 = [STAs(:)/norm(STAs(:)), PCs(:,1), PCs(:,2), PCs(:,3), PCs(:,4)];
    new_vec = [tmp1(:,str2num(ch1)), tmp1(:,str2num(ch2))]; % normalising the vectors

else
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
        return_basis_vec = PCs(:,choice-1);
        return_basis_vec = reshape(return_basis_vec,[nrandnums_perchannel*3 maxT]);
    end
end

 
%% A different analysis starts here
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
    
    figure_label = ['Time Course: STA';'Time Course: PC1';'Time Course: PC2'];
    vectors = [STAs(:), PCs(:,1), PCs(:,2)];
    for j = 1:size(vectors,2)
        vec_temp = reshape(vectors(:,j),[nrandnums_perchannel*3 maxT]);
        for i = 1:subunits
            temp_plot(1,:,i) = vec_temp(i,:); % R
            temp_plot(2,:,i) = vec_temp(subunits+i,:); % G
            temp_plot(3,:,i) = vec_temp(2*subunits+i,:); % B
            tmp = squeeze(temp_plot(:,:,i));
            LMS_temp_plot(:,:,i) = M*tmp;
            
            % plotting RGB signal variation across frames
            h1 = figure(plot_counter);
            set(h1,'name','RGB');subplot(size(vectors,2),subunits,(j-1)*subunits+i);
            plot(temp_plot(1,:,i),'r','LineWidth',2); hold on;
            plot(temp_plot(2,:,i),'g','LineWidth',2); hold on;
            plot(temp_plot(3,:,i),'LineWidth',2); hold on;
            xlabel('frames'),ylabel('Intensity'); title(figure_label(j,:)); hold off;
            
            % plotting LMS signal variation across frames
            h2 = figure(plot_counter+1);
            set(h2,'name','LMS'); subplot(size(vectors,2),subunits,(j-1)*subunits+i);
            plot(LMS_temp_plot(1,:,i),'r','LineWidth',2); hold on;
            plot(LMS_temp_plot(2,:,i),'g','LineWidth',2); hold on;
            plot(LMS_temp_plot(3,:,i),'LineWidth',2); hold on;
            xlabel('frames'),ylabel('Intensity');title(figure_label(j,:)); hold on;
        end
        
        clear vec_temp tmp
    end
end
hold off;
plot_counter = plot_counter + 2;

% Displaying the basis vectors
FigHandle = figure(plot_counter);
set(FigHandle, 'Position', [50, 250, 1449, 350]);
if flag == 0 
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
else
        
    vec1 = new_vec(:,1); vec2 = new_vec(:,2);
    vec1 = reshape(vec1,[nrandnums_perchannel*3 maxT]);
    vec2 = reshape(vec2,[nrandnums_perchannel*3 maxT]);
    for i = 1:maxT % evaluates for each frame
        normfactor1 = 0.5/((max(abs(vec1(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
        tmp1 = vec1(:,i);
        tmp1 = expand_vector(tmp1,nrandnums_perchannel,mask,1);
        tmp1 = normfactor1*(tmp1)+muvect; % This makes the values fall back within a range of 0 and 1.
        tmp1 = reshape(tmp1,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
        
        normfactor2 = 0.5/((max(abs(vec2(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5
        tmp2 = vec2(:,i);
        tmp2 = expand_vector(tmp2,nrandnums_perchannel,mask,1);
        tmp2 = normfactor2*(tmp2)+muvect; % This makes the values fall back within a range of 0 and 1.
        tmp2 = reshape(tmp2,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
        
        subplot(2,maxT,i);
        image(abs(tmp1)); % for viewing the image
        set(gca,'XTick',[],'YTick',[]); axis square;
        if (i == 1)
            if (str2num(ch1) == 1)
                ylabel('STA');
            elseif (str2num(ch1) == 2)
                ylabel('PC1');
            elseif (str2num(ch1) == 3)
                ylabel('PC2');
            elseif (str2num(ch1) == 4)
                ylabel('PC3');
            elseif (str2num(ch1) == 5)
                ylabel('PCi');
            end
        end
        
        subplot(2,maxT,maxT+i);
        image(abs(tmp2)); % for viewing the image
        set(gca,'XTick',[],'YTick',[]); axis square;
        
        if (i == 1)
             if (str2num(ch2) == 1)
                 ylabel('STA');
             elseif (str2num(ch2) == 2)
                 ylabel('PC1');
             elseif (str2num(ch2) == 3)
                 ylabel('PC2');
             elseif (str2num(ch2) == 4)
                 ylabel('PC3');
             elseif (str2num(ch2) == 5)
                 ylabel('PCi');
             end
        end
    end
end
    
plot_counter = plot_counter + 1;

end

