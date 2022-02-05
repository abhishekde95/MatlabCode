% Post-rebutall script
% Author -Abhishek De
% Date - 02/22

close all; clearvars;

%% 1. Relationship between Isoresponse NLI and target firing rates

if ~exist('plot_counter')
    plot_counter = 1;
end

% Loading all the files 
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end


load conewts_svd.mat
load vals.mat
load vals.mat
load S1RGB_svd.mat
load S2RGB_svd.mat
load angulardifferences_RGB.mat
anglebwvectors = angulardifference_RGB;
S1RGB = S1RGB_svd;
S2RGb = S2RGB_svd;
% SpatiallyOpponent = sum(sign(S1RGB).*sign(S2RGB),1)<3;
SpatiallyOpponent = anglebwvectors'>90;

thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Considering only the spatially opponent subunits
LUMidx = LUMidx(SpatiallyOpponent(LUMidx));
DOidx = DOidx(SpatiallyOpponent(DOidx));
hardtoclassifyidx = hardtoclassifyidx(SpatiallyOpponent(hardtoclassifyidx));

% Checking how target firing rates of the DO, simple and unclassified cells
load baselineFRstats.mat % Baseline FR rates
load TFR.mat % Target firing rates


% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat


% For storing the Isoresponse NLI
Isoresponse_NLI = [];


for ii = 1:numel(RSSE_linearmodel) 
   
    % Isoresponse NLI
    Isoresponse_NLI = [Isoresponse_NLI; log10(median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii}))];
end

%%
idx = [LUMidx, DOidx, hardtoclassifyidx];
[r, p] = corr(Isoresponse_NLI(idx), TFR(1,idx)', 'type', 'Spearman');
