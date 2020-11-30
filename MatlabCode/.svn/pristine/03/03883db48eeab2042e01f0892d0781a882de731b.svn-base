%% This program runs a series of files to fit a psychometric function to data.
clear all
close all

%% User-defined Variables

numdatafiles=2;
datafile1= 'SK090211';
datafile2= 'SK080211';
%datafile1 = 'JPW110214';
%datafile2 = 'JPW110215';
num_resamples=250;
tg=[1 1 1 1 1 1 1];


%% Auto-defined Variables
datafile1=([datafile1 '.xml']);
if numdatafiles==2
    datafile2=([datafile2 '.xml']);
    datafiles=[datafile1; datafile2];
else
    datafiles=datafile1;
end


%% Loop process to analyze multiple datafiles
for q=1:size(datafiles,1)
    
    % Selects info from an .xml file and converts into a .mat file
    [datafile_eval1(:,:,q)] = xml2mat(datafiles(q,:));
    
    %Reorganizes rawdata and plots it
    [df(q)] = orgrawdata(datafile_eval1(:,:,q),q);
    
    %Fits curves to both psychometric and chronometric data
    [th_fit(q,:), err(:,q)] = fitwithcurves(df(q),tg(1:4),q)
   
    % Resamples original data, fits curves, and plots range for each varable
    [th_fit_rs(:,:,q),err_rs(:,q),rs(1,:,q)] = resample(datafile_eval1(:,:,q),tg(1:4),num_resamples);
    
    % Sanity checks for recovering set variables
    n = 5000;
    
    % Sanity check for recovering different k's
    [th_fit_san,err_san] = sanitycheck(df,tg,n)
    %[th_fit_sant,err_sant] = sanitycheck_t2(df,tg,n)
    
    % Sanity check for recovering different A's
    [th_fit_san2,err_san2] = sanitycheck2(df,tg,n)
        
    % Sanity check for recovering different B's
    [th_fit_san3,err_san3] = sanitycheck3(df,tg,n)
    
end

%% Fit datasets together by high and low coherences
% One A and Tnd, no B, but different k's for H and L coherences
[th_fit1,err1] = fitwithcurves1(df,tg(1:4))

% One A, Tnd, and B, but different k's for H and L coherences
[th_fit2,err2] = fitwithcurves2(df,tg(1:5))
%[th_fit2_t,err2_t] = fitwithcurves2_t(df,tg(1:5))

% One A, Tnd, and k, but different B's for H and L coherences
[th_fit3,err3] = fitwithcurves3(df,tg(1:5))

% One k, Tnd, and B, but different A's for H and L coherences
[th_fit4,err4] = fitwithcurves4(df,tg(1:5))

% One k, Tnd, B, and A for H and L coherences
[th_fit5,err5] = fitwithcurves5(df,tg(1:4))

% One k and Tnd, no B, and different A's for H and L coherences
[th_fit6,err6] = fitwithcurves6(df,tg(1:4))

% Two k's, two A's, one Tnd and one B for H and L coherences
[th_fit7,err7] = fitwithcurves7(df,tg(1:6))

% Resample 7 using rs data
[kHdist7,kLdist7,AHdist7,ALdist7] = resample_th_fit7(rs,tg(1:6),num_resamples)
        
% Two k's, two A's, one Tnd and no B
[th_fit8,err8] = fitwithcurves8(df,tg(1:5))

% Two k's, two A's, two B's, and one Tnd
[th_fit9,err9] = fitwithcurves9(df,tg(1:7))

% Two k's, two B's, one A, and one Tnd
[th_fit10,err10] = fitwithcurves10(df,tg(1:6))

