% Author - Abhishek De, 11/16
% Script for understanding linear/nonlinear spatial integration 
% The idea is to look at if spatial nonlinear combination can give rise to
% non-planar isoresponse surfaces in the gun/cone space
clearvars; close all;
% Create an STA of a hypothetical neuron, RF size L x W pixels
L = 20; W = 20;
mask = zeros(L,W);
mask(round(L/4):round(3*L/4),round(W/4):round(2*W/4)) = 1; % Subunit 1
mask(round(L/4):round(3*L/4),round(2*W/4)+1:round(3*W/4)) = 2; % Subunit 2, rest all is background
bg = find(mask == 0);
S1 = find(mask == 1);
S2 = find(mask == 2);
RF = zeros(L*W,3);
mode = 0;
if mode == 0
    load RGBS1.mat;
    load RGBS2.mat;
elseif mode == 1
    RGBS1  = randn(3,1); RGBS2  = randn(3,1);
elseif mode == 2
    RGBS1  = 0.3*ones(3,1); RGBS2  = -0.3*ones(3,1);
end
RF(S1,1) = RGBS1(1); RF(S1,2) = RGBS1(2); RF(S1,3) = RGBS1(3);
RF(S2,1) = RGBS2(1); RF(S2,2) = RGBS2(2); RF(S2,3) = RGBS2(3);

% Create a gabor template
theta = pi/2; sigma = 0.4; ggamma = 1; sf = 0.5;
[x,y] = meshgrid(linspace(-1,1,L), linspace(-1,1,W));
X = x*cos(-theta) + y*sin(-theta);
Y =-x*sin(-theta) + y*cos(-theta);
expterm = exp(-(X.^2 + ggamma^2*Y.^2)/2/sigma^2);
phasecount = 0;

% Start running the simulation 
spatial_summation = 4; % 1 - linear; 2 - squaring and adding; 3 - hyperbolic; 4 - squaring and subtracting
for phase = 0 % As of now just implement for 1 particular phase
    phasecount = phasecount + 1;
    costerm = cos(2*pi*Y*sf + phase);
    template = costerm.*expterm;
     
    % Start with 6 initial directions and run it through the adaptive closed loop algo
    dirs = [1 0 0; 0 1 0; 0 0 1; -1 0 0;0 -1 0; 0 0 -1];
    N = size(dirs,1);
    contrast = repmat(0.5,[N 1]);
    reversals = zeros([N 1]);
    active = ones([N 1]);
    done = zeros([N 1]);
    parents = [];
    dirs_to_probe = 200;
    target_FR = 40;
    nreversals = 20;
    scale_contrast = repmat(0.25,[N 1]);
    end_FR = zeros([N 1]);
    oog = zeros([N 1]);
    gamut_edge = 1;
    scale_fact = 0.75;
    contrast_traj = cell(N,1);
    FR_traj = cell(N,1);
    Out_mat = [];
    driveS1_mat = [];
    driveS2_mat = [];
    
    while numel(contrast)<dirs_to_probe
        ind = find(done==0);
        if isempty(ind) % When all the directions in the queue have been probed
            if isempty(parents)
                active(:) = 0;
                new_ind = numel(active)+1;
                parents = [0 0 0; 0 0 0; 0 0 0;0 0 0;0 0 0;0 0 0; 1 2 3; 1 2 6;4 2 3; 4 2 6; 4 5 3; 4 5 6; 1 5 3; 1 5 6];
                reversals(new_ind:new_ind+7) = zeros([8 1]);
                done(new_ind:new_ind+7) = zeros([8 1]);
                active(new_ind:new_ind+7) = ones([8 1]);
                contrast(new_ind:new_ind+7) = repmat(0.5,[8 1]);
                scale_contrast(new_ind:new_ind+7) = repmat(0.25,[8 1]);
                end_FR(new_ind:new_ind+7) = zeros([8 1]);
                oog(new_ind:new_ind+7) = zeros([8 1]);
                new_dir1 = contrast(1)*dirs(1,:)+ contrast(2)*dirs(2,:)+ contrast(3)*dirs(3,:);
                new_dir2 = contrast(1)*dirs(1,:)+ contrast(2)*dirs(2,:)+ contrast(6)*dirs(6,:);
                new_dir3 = contrast(4)*dirs(4,:)+ contrast(2)*dirs(2,:)+ contrast(3)*dirs(3,:);
                new_dir4 = contrast(4)*dirs(4,:)+ contrast(2)*dirs(2,:)+ contrast(6)*dirs(6,:);
                new_dir5 = contrast(4)*dirs(4,:)+ contrast(5)*dirs(5,:)+ contrast(3)*dirs(3,:);
                new_dir6 = contrast(4)*dirs(4,:)+ contrast(5)*dirs(5,:)+ contrast(6)*dirs(6,:);
                new_dir7 = contrast(1)*dirs(1,:)+ contrast(5)*dirs(5,:)+ contrast(3)*dirs(3,:);
                new_dir8 = contrast(1)*dirs(1,:)+ contrast(5)*dirs(5,:)+ contrast(6)*dirs(6,:);
                dirs(new_ind,:) = new_dir1/norm(new_dir1); dirs(new_ind+1,:) = new_dir2/norm(new_dir2);
                dirs(new_ind+2,:) = new_dir3/norm(new_dir3); dirs(new_ind+3,:) = new_dir4/norm(new_dir4);
                dirs(new_ind+4,:) = new_dir5/norm(new_dir5); dirs(new_ind+5,:) = new_dir6/norm(new_dir6);
                dirs(new_ind+6,:) = new_dir7/norm(new_dir7); dirs(new_ind+7,:) = new_dir8/norm(new_dir8);
                contrast_traj{new_ind} = []; contrast_traj{new_ind+1} = [];contrast_traj{new_ind+2} = []; contrast_traj{new_ind+3} = [];
                contrast_traj{new_ind+4} = []; contrast_traj{new_ind+5} = [];contrast_traj{new_ind+6} = []; contrast_traj{new_ind+7} = [];
                FR_traj{new_ind} = []; FR_traj{new_ind+1} = [];FR_traj{new_ind+2} = []; FR_traj{new_ind+3} = [];
                FR_traj{new_ind+4} = []; FR_traj{new_ind+5} = [];FR_traj{new_ind+6} = []; FR_traj{new_ind+7} = [];
                disp('Creating the fourth direction');
            else
                active_ind = find(active == 1);
                active(:) = 0;
                clear bb;
                for bb = 1:numel(active_ind)
                    jj = active_ind(bb);
                    grandparents = [parents(jj,:) parents(jj,1)];
                    for kk = 1:3
                        new_ind = numel(active) + 1;
                        parents(new_ind,:) = [jj grandparents(kk) grandparents(kk+1)];
                        reversals(new_ind) = 0;
                        done(new_ind) = 0;
                        active(new_ind) = 1;
                        contrast(new_ind) = 0.5;
                        scale_contrast(new_ind) = 0.25;
                        end_FR(new_ind) = 0;
                        oog(new_ind) = 0;
                        new_dir = 0.5*dirs(jj,:)+ 0.5*dirs(grandparents(kk),:)+ 0.5*dirs(grandparents(kk+1),:);
                        dirs(new_ind,:) = new_dir/norm(new_dir);
                        contrast_traj{new_ind} = [];
                        FR_traj{new_ind} = [];
                    end
                end
                disp('Creating new directions');
            end
            ind = find(done==0);
        end
        for aa = 1:numel(ind)
            ii = ind(aa);
            Inp = [dirs(ii,1)*template(:) dirs(ii,2)*template(:) dirs(ii,3)*template(:)];
            Inp = contrast(ii)*Inp;
            count = 0;
            prev_out = [];
            while reversals(ii)<nreversals
                [Out,driveS1,driveS2,net_drive] = get_model_spatial_response(Inp*contrast(ii),RF,S1,S2,spatial_summation);
                Out_mat = [Out_mat; Out];
                driveS1_mat = [driveS1_mat; driveS1];
                driveS2_mat = [driveS2_mat; driveS2];
                count = count + 1;
                if count == 1
                    prev_out = Out;
                else
                    if (Out-target_FR)*(prev_out-target_FR)<=0
                        reversals(ii) = reversals(ii)+1;
                        scale_contrast(ii) = scale_contrast(ii)*scale_fact;
                    end
                end
                if Out<=target_FR
                    contrast(ii) = contrast(ii)*(1+scale_contrast(ii));
                else
                    contrast(ii) = contrast(ii)*(1-scale_contrast(ii));
                end
                FR_traj{ii} = [FR_traj{ii};Out];
                contrast_traj{ii}= [contrast_traj{ii};contrast(ii)];
                tmp_dir = dirs(ii,:)*contrast(ii);
                if (contrast(ii)>gamut_edge)
                    oog(ii) = 1;
                    break;
                end
            end
            done(ii) = 1;
            end_FR(ii) = FR_traj{ii}(end);
        end
    end
    tmp = repmat(contrast,[1 3]);
    new_dirs = dirs.*tmp;
end
disp('Done running the closed loop algo');

% More Analysis
tmp_idx = find(oog==0); % Only the the inside gamut points
tmp = repmat(contrast,[1 3]);
new_dirs = dirs.*tmp;
figure(3),
subplot(231),image(XW2RGBFormat(RF,L,W)+0.5); title('STA of the neuron');
subplot(232),scatter3(new_dirs(tmp_idx,1),new_dirs(tmp_idx,2),new_dirs(tmp_idx,3),10,'filled'); title('Isoresponse contour'); view([30,20]);
if spatial_summation
    subplot(232),hold on; scatter3(-new_dirs(tmp_idx,1),-new_dirs(tmp_idx,2),-new_dirs(tmp_idx,3),10,'filled'); title('Isoresponse contour'); view([30,20]); hold off;
end
subplot(233),scatter3(dirs(:,1),dirs(:,2),dirs(:,3),'r','filled'), title('Unit vectors'),view([30,20]);
subplot(234),imagesc(template);
subplot(235); hist(contrast); xlabel('Contrast'),ylabel('Frequency');
subplot(236), hist(end_FR); xlabel('Firing rate'),ylabel('Frequency');


