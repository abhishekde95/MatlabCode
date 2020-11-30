% This script is to understand how spatial summation happens for cells
% which have a spatial structure in STA and PC1. I am doing this to mostly
% come up with predictions
% Author - Abhishek De, 10/17
close all; clearvars;
L = 100; % pixels
STA = zeros(L,L,3); % 3 is for L,M,S channels
PC = zeros(L,L,3);
% Create a cell that has cone-opponent STA and cone-non-opponent PC1
% subunit 1
RGB = 0.1*[0.25 0.25 -0.5];%rand(3,1);
RGB1 = [0.5 0.5 0.5];
STA(20:70,20:50,1) = RGB(1); STA(20:70,20:50,2) = RGB(2); STA(20:70,20:50,3) = RGB(3); % M goes along with S
subunit1 = STA; 
PC(20:70,20:50,1) = RGB1(1); PC(20:70,20:50,2) = RGB1(2); PC(20:70,20:50,3) = RGB1(3);
% subunit 2
STA(20:70,51:80,1) = RGB(1); STA(20:70,51:80,2) = RGB(2); STA(20:70,51:80,3) = RGB(3); % M goes along with S
subunit2 = STA - subunit1; %subunit2 = subunit2./norm(subunit2(:));
%subunit1 = subunit1./norm(subunit1(:));
PC(20:70,51:80,1) = -RGB1(1); PC(20:70,51:80,2) = -RGB1(2); PC(20:70,51:80,3) = -RGB1(3);
% STA = STA./norm(STA(:));
% PC = PC./norm(PC(:));
figure(1),subplot(221),image((0.5*STA./(max(abs(STA(:))))) + 0.5); set(gca,'XTick',[],'YTick',[]); title('STA');
subplot(222),image((0.5*PC./(max(abs(PC(:))))) + 0.5); set(gca,'XTick',[],'YTick',[]); title('PC');
subplot(223),image((0.5*subunit1./(max(abs(subunit1(:))))) + 0.5); set(gca,'XTick',[],'YTick',[]); title('Basis Vec 1');
subplot(224),image((0.5*subunit2./(max(abs(subunit2(:))))) + 0.5); set(gca,'XTick',[],'YTick',[]); title('Basis Vec 2');

% Step 1 is done, have created the filters, next step is to carry out the
% isoresponse closed loop measuremets on this hypothetical cell
% Start with 5 initial directions and run it through the adaptive closed
% loop algo, pretty similar to what I do in my experiments
dirs = [-sqrt(1/2) sqrt(1/2); 0 1; sqrt(1/2) sqrt(1/2); 1 0; sqrt(1/2) -sqrt(1/2)];
N = size(dirs,1);
contrast = repmat(0.5,[N 1]);
reversals = zeros([N 1]);
active = ones([N 1]);
done = zeros([N 1]);
parents = [];
dirs_to_probe = 50;
target_FR = 30;
nreversals = 20;
scale_contrast = repmat(0.25,[N 1]);
end_FR = zeros([N 1]);
oog = zeros([N 1]);
gamut_edge = 2; % Weber contrast
scale_fact = 0.75;
contrast_traj = cell(N,1);
FR_traj = cell(N,1);
Out_mat = [];

while numel(contrast)<dirs_to_probe
    numel(contrast)
    ind = find(done==0);
    if isempty(ind) % When all the directions in the queue have been probed
        if isempty(parents)
            active(:) = 0;
            new_ind = numel(active)+1;
            parents = [0 0; 0 0; 0 0; 0 0; 0 0; 1 2; 2 3; 3 4; 4 5];
            reversals(new_ind:new_ind+3) = zeros([4 1]);
            done(new_ind:new_ind+3) = zeros([4 1]);
            active(new_ind:new_ind+3) = ones([4 1]);
            contrast(new_ind:new_ind+3) = repmat(0.5,[4 1]);
            scale_contrast(new_ind:new_ind+3) = repmat(0.25,[4 1]);
            end_FR(new_ind:new_ind+3) = zeros([4 1]);
            oog(new_ind:new_ind+3) = zeros([4 1]);
            new_dir1 = contrast(1)*dirs(1,:)+ contrast(2)*dirs(2,:);
            new_dir2 = contrast(2)*dirs(2,:)+ contrast(3)*dirs(3,:);
            new_dir3 = contrast(3)*dirs(3,:)+ contrast(4)*dirs(4,:);
            new_dir4 = contrast(4)*dirs(4,:)+ contrast(5)*dirs(5,:);
            dirs(new_ind,:) = new_dir1/norm(new_dir1); dirs(new_ind+1,:) = new_dir2/norm(new_dir2);
            dirs(new_ind+2,:) = new_dir3/norm(new_dir3); dirs(new_ind+3,:) = new_dir4/norm(new_dir4);
            contrast_traj{new_ind} = []; contrast_traj{new_ind+1} = [];
            contrast_traj{new_ind+2} = []; contrast_traj{new_ind+3} = [];
            FR_traj{new_ind} = []; FR_traj{new_ind+1} = [];
            FR_traj{new_ind+2} = []; FR_traj{new_ind+3} = [];
            disp('Creating the next set of directions');
        else
            active_ind = find(active == 1);
            active(:) = 0;
            clear bb;
            for bb = 1:numel(active_ind)
                jj = active_ind(bb);
                grandparents = parents(jj,:);
                for kk = 1:2
                    new_ind = numel(active) + 1;
                    parents(new_ind,:) = [jj grandparents(kk)];
                    reversals(new_ind) = 0;
                    done(new_ind) = 0;
                    active(new_ind) = 1;
                    contrast(new_ind) = 0.5;
                    scale_contrast(new_ind) = 0.25;
                    end_FR(new_ind) = 0;
                    oog(new_ind) = 0;
                    new_dir = 0.5*dirs(jj,:)+ 0.5*dirs(grandparents(kk),:);
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
        Inp = dirs(ii,1)*subunit1 + dirs(ii,2)*subunit2;
        count = 0;
        prev_out = [];
        while reversals(ii)<nreversals
            response = getresponse(contrast(ii)*Inp,STA,PC,subunit1,subunit2);
            Out_mat = [Out_mat; response];
            count = count + 1;
            if count == 1
                prev_out = response;
            else
                if (response-target_FR)*(prev_out-target_FR)<=0
                    reversals(ii) = reversals(ii)+1;
                    scale_contrast(ii) = scale_contrast(ii)*scale_fact;
                end
            end
            if response<=target_FR
                contrast(ii) = contrast(ii)*(1+scale_contrast(ii));
            else
                contrast(ii) = contrast(ii)*(1-scale_contrast(ii));
            end
            FR_traj{ii} = [FR_traj{ii}; response];
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
endpts = repmat(contrast,[1 2]).*dirs;
[THETA,RHO] = cart2pol(endpts(:,1),endpts(:,2));

% Trying out a new thing just to confirm if I have the simulation correctly
STAcoeffs = linspace(-2,2,40);
PCcoeffs = linspace(-2,2,40);
subunit1coeffs = STAcoeffs;
subunit2coeffs = STAcoeffs;
[STAcoeffs, PCcoeffs] = meshgrid(STAcoeffs,PCcoeffs);
[subunit1coeffs,subunit2coeffs] = meshgrid(subunit1coeffs,subunit2coeffs);
respmatSTAPC = zeros(size(STAcoeffs));
respmatsubunit = zeros(size(STAcoeffs));
projontoSTA = [];
projontoPC = [];
projontosubunit1 = [];
projontosubunit2 = [];
for ii = 1:size(STAcoeffs,1)
    for jj = 1:size(STAcoeffs,2)
        Inp1 = STAcoeffs(ii,jj)*STA + PCcoeffs(ii,jj)*PC;
        respmatSTAPC(ii,jj) = getresponse(Inp1,STA,PC,subunit1,subunit2);
        Inp2 = subunit1coeffs(ii,jj)*subunit1 + subunit2coeffs(ii,jj)*subunit2;
        respmatsubunit(ii,jj) = getresponse(Inp2,STA,PC,subunit1,subunit2);
        projontosubunit1 = [projontosubunit1; Inp1(:)'*subunit1(:)];
        projontosubunit2 = [projontosubunit2; Inp1(:)'*subunit2(:)];
        projontoSTA = [projontoSTA; Inp2(:)'*STA(:)];
        projontoPC = [projontoPC; Inp2(:)'*PC(:)];
    end
end

initpts = repmat([0 0],[sum(oog) 1]);
oogidx = find(oog==1);
figure(2),subplot(231), hist(contrast); ylabel('Frequency'); xlabel('Final Contrast'); 
subplot(232), plot(endpts(~oog,1),endpts(~oog,2),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'PickableParts','none','MarkerEdgeColor',[1 0 0]); hold on;
for kk = 1:size(initpts,1)
    plot([initpts(kk,1) endpts(oogidx(kk),1)],[initpts(kk,2) endpts(oogidx(kk),2)],'g'); % Drawing the out of gamut points
end
set(gca,'Xlim',[-1.5 1.5],'Ylim',[-1.5 1.5]); plot([-1 1],[0 0],'-k','Linewidth',2); plot([0 0],[-1 1],'-k','Linewidth',2); xlabel('Basis Vec 1'); ylabel('Basis Vec 2');
title('Isoresponse contour');hold off;
subplot(234),hist(end_FR); ylabel('Frequency'), xlabel('Final FR');
subplot(233),contour3(subunit1coeffs,subunit2coeffs,respmatsubunit,'LineWidth',2); xlabel('subunit1'), ylabel('subunit2'); zlabel('FR'); title('subunit combination'); set(gca,'XTick',[],'YTick',[]);
subplot(235),contour3(STAcoeffs,PCcoeffs,respmatSTAPC,'LineWidth',2); xlabel('STA'), ylabel('PC'); zlabel('FR'); title('STA vs PC combination'); set(gca,'XTick',[],'YTick',[]);

% Trying to see how the cone signal combination takes place 
theta = pi/2; sigma = 0.4; ggamma = 1; sf = 1.0;
[x,y] = meshgrid(linspace(-1,1,L), linspace(-1,1,L));
X = x*cos(-theta) + y*sin(-theta);
Y =-x*sin(-theta) + y*cos(-theta);
expterm = exp(-(X.^2 + ggamma^2*Y.^2)/2/sigma^2);
phasecount = 0;

% Start running the simulation 
cone_signal_combination = 1; % 1- linear; 2 - ellipsoid, 3 - hyperboloid
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
    dirs_to_probe = 100;
    nreversals = 25;
    scale_contrast = repmat(0.25,[N 1]);
    end_FR = zeros([N 1]);
    target_FR = 30;
    oog = zeros([N 1]);
    gamut_edge = 5;
    scale_fact = 0.75;
    contrast_traj = cell(N,1);
    FR_traj = cell(N,1);
    Out_mat = [];
    
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
            Inp = cat(3,dirs(ii,1)*template,dirs(ii,2)*template,dirs(ii,3)*template);
            
            count = 0;
            prev_out = [];
            while reversals(ii)<nreversals
%                 contrast(ii)
                Out = getresponse(contrast(ii)*Inp,STA,PC,subunit1,subunit2);
                Out_mat = [Out_mat; Out];
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
oogidx = find(oog==1);
tmp = repmat(contrast,[1 3]);
new_dirs = dirs.*tmp;
figure(2), subplot(236), plot3(new_dirs(~oog,1),new_dirs(~oog,2),new_dirs(~oog,3),'o','LineWidth',1.0,'MarkerFaceColor',[1 0 0],'MarkerSize',6); hold on;
for kk = 1:numel(oogidx)
    plot3([0 new_dirs(oogidx(kk),1)],[0 new_dirs(oogidx(kk),2)],[0 new_dirs(oogidx(kk),3)],'g'); % Drawing the out of gamut points
end
xlabel('L'), ylabel('M'),zlabel('S'), title('Cone signal combination');
plot3([-3 3],[0 0],[0 0],'-k','Linewidth',2); % x- axis 
plot3([0 0],[-3 3],[0 0],'-k','Linewidth',2); % y- axis
plot3([0 0],[0 0],[-1 1],'-k','Linewidth',2); % z- axis
hold off;
