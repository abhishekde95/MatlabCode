%This seems like an old part of lmtf_generate_module_data, or maybe
%something to do with later model fitting - probably replaced by
%lmtf_global_to_local_model
function directModelParams = fitModelFromData(subj)
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
f_q = sprintf('CALL postPropixxFilenames(''%s'');', subj);
filenames = fetch(conn, f_q);
[data, ~] = iterateAndPlotFiles_modularPlusDB(filenames, 0);
LB = MakeParamList([bounds.zeta(1) bounds.n(1) bounds.delta_n(1) bounds.logtau(1) bounds.logkappa(1) bounds.zeta(1) bounds.n(1) bounds.delta_n(1) bounds.logtau(1) bounds.logkappa(1) bounds.theta(1)],...
                   [bounds.a1(1) bounds.a1(1) bounds.a1(1) bounds.a1(1) bounds.a1(1) bounds.a1(1) bounds.a1(1) bounds.a1(1)]);
UB = MakeParamList([bounds.zeta(2) bounds.n(2) bounds.delta_n(2) bounds.logtau(2) bounds.logkappa(2) bounds.zeta(2) bounds.n(2) bounds.delta_n(2) bounds.logtau(2) bounds.logkappa(2) bounds.theta(2)],...
                   [bounds.a1(2) bounds.a1(2) bounds.a1(2) bounds.a1(2) bounds.a1(2) bounds.a1(2) bounds.a1(2) bounds.a1(2)]);
%initial guess goes here
[fpar,~] = fmincon(@(params) tf_fiterr3(params,data,4),initialguess,[],[],[],[],LB,UB,[],options); % mode 4 fitting (double tilted rampy trough)
directModelParams = fpar;
end