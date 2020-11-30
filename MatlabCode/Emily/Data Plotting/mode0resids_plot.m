%plots 11+2n model 

function mode0resids_plot(subjID)
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
f_req = sprintf('CALL postPropixxFilenames(''%s'')', subjID);
flist = fetch(conn, f_req);
fulldata = getLMTFrawdata(flist);
try
    struct = load('C:\Documents and Settings\emily\Desktop\MatlabCode\trunk\Zack\IsoSamp\private\data\LMTF.mat', subjID);
catch
    struct = load('C:\Users\emily.gelfand\Desktop\MatlabCode\trunk\Zack\IsoSamp\private\data\LMTF.mat', subjID);
end
model = getfield(struct, subjID);
global_model = model.legacy.mode0models;
model_eccs = model.eccs;
eccs = unique(fulldata(:,[5 6]), 'rows');
figure; hold on; title([subjID ' 11+2n']);
all_resids = []; lengths = [];
for i = 1:length(eccs)
    Lecc = fulldata(:,5)==eccs(i,1) & fulldata(:,6)==eccs(i,2);
    ecc_Data = fulldata(Lecc, :);
    Loog = ecc_Data(:,4);
    L = ecc_Data(~Loog,1);
    M = ecc_Data(~Loog,2);
    TF = ecc_Data(~Loog,3);
    model_col_idx = find(model_eccs(:,1) == eccs(i,1) & model_eccs(:,2) == eccs(i,2));
    local_model = global_model(:,model_col_idx);
    predr = LMTF_thresh_from_model(L, M, TF, local_model);
    realr = sqrt(L.^2 + M.^2);
    resids = log(realr) - log(predr);
    all_resids = [all_resids resids'];
    lengths = [lengths length(resids)];
    labels{i} = num2str(model_eccs(model_col_idx,:));
end
g = [];
for j = 1:length(lengths)
   if i==1
       g = zeros(1,lengths(j));
   else
       g = [g (j-1)*ones(1,lengths(j))];
   end
end
p_val = kruskalwallis(all_resids, g, 'on');
disp(p_val);
boxplot(all_resids, g, 'labels', labels);