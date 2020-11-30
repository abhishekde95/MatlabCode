% This is a population analysis script is adopted from Greg's Grating 
% Population Analysis.  Reading datafile names from a .txt file, it will 
% organize the data into a giant cell array.  Once in this format,
% poulation analyses can be performed.

clear all
close all


%% Load all datafiles

if isdir('Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles')
    [fnames, spikeIdx] = fnamesFromTxt2('Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/GLMPnexfilelist.txt');
    startingPath = 'Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles';
elseif isdir('N:\NexFiles\Patrick\')
    [fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\Patrick\GLMPnexfilelist.txt');
    startingPath = 'N:\NexFiles\Patrick\';
else
    disp('Nothing in path...')
    keyboard
end

data = cell(size(fnames,1),1);

for a = 1:size(fnames,1)
    stro = {};
    for i = 1:size(fnames{a},2)
        disp(['Retrieving datafile ' char(fnames{a}(i)) '...'])
        tmpstro = nex2stro(findfile(char(fnames{a}(i)),startingPath));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end

    disp(['Unpacking datafile ' char(stro.sum.fileName(end-13:end-4)) '...']);
    [trial par plat] = OrganizeRawGLMPData(stro);
    data{a} = {trial par plat};
    
end


%% Set up a few variables

fig = 900;


%% Begin Individual Analyses
% 1. Anova1

for d = 1:size(data,1)
    
    datafile = data{d}{1}.datafile;
    disp('Anova1')
    disp(['Now processing ' datafile])
    
    
    if d == 1
        analyzePlats = [];
        platNames = [];
        platNum = [];
    end
    
    data{d}{3} = GLMP_AOV(data{d}{3});


    for p = 1:numel(data{d}{3})
        
        % Do not include adaptive files
        if cell2mat(data{d}{3}{p}.allProp(end,2)) == 1
            data{d}{3}{p}.Analyze = 0;
        end
        
        analyzePlats = cat(1,analyzePlats,data{d}{3}{p}.Analyze);
        platNames = cat(1,platNames,data{d}{3}{p}.datafile);
        platNum = cat(1,platNum,p);
        
    end
    
    disp('Should exclude the following datafiles:')
    
    platNames(analyzePlats==0,:)
    platNum(analyzePlats==0)
    
end
    
  
%% 1-D Fit to All Data
% 
% bestFitAngs1D_AllData = [];
% 
% for d = 1:size(data,1)
%     
%     datafile = data{d}{1}.datafile;
%     disp('1-D Fit To All Data')
%     disp(['Now processing ' datafile])
%     
%     data{d}{3} = OneDFit_AllData(data{d}{3});
%     
%     for p = 1:numel(data{d}{3})
%         
%         if data{d}{3}{p}.Analyze == 1
%             
%             tempAng = data{d}{3}{p}.OneDFitAllData.Best.RotationInDeg;
%             bestFitAngs1D_AllData = cat(1,bestFitAngs1D_AllData,tempAng);
%             
%         end
%         
%     end
% 
% end
% 
% figure(fig); clf; hold on; grid on; fig=fig+1;
% plotTitle = 'Population 1-D Fit to All Data';
% set(gcf,'Name',plotTitle,'NumberTitle','off')
% hist(bestFitAngs1D_AllData,32)
% xlabel('Best Fit Rotation (deg)')
% ylabel('Number of Cells')
% 

     

%% 2-D Fit to All Data

bestFitAngs2D_AllData = [];

for d = 1:size(data,1)
    
    datafile = data{d}{1}.datafile;
    disp('2-D Fit To All Data')
    disp(['Now processing ' datafile])
    
%     try
%         data{d}{3} = TwoDFit_AllData(data{d}{3},1);
%     catch
%         data{d}{3} = TwoDFit_AllData2(data{d}{3},1);
%     end
    data{d}{3} = TwoDFit_AllData2(data{d}{3},1);
            
        
    for p = 1:numel(data{d}{3})
        
        if data{d}{3}{p}.Analyze == 1
            
            tempAng = data{d}{3}{p}.TwoDFitAllData.Best.RotationInDeg;
            bestFitAngs2D_AllData = cat(1,bestFitAngs2D_AllData,tempAng);
            
        end
        
    end
    
end

%A bit of formatting
bestFitAngs2D_AllData(bestFitAngs2D_AllData<0) = bestFitAngs2D_AllData(bestFitAngs2D_AllData<0)+360;
bestFitAngs2D_AllData(bestFitAngs2D_AllData>90 & bestFitAngs2D_AllData<=180)=...
    bestFitAngs2D_AllData(bestFitAngs2D_AllData>90 & bestFitAngs2D_AllData<=180)-90;
bestFitAngs2D_AllData(bestFitAngs2D_AllData>180 & bestFitAngs2D_AllData<=270)...
    =bestFitAngs2D_AllData(bestFitAngs2D_AllData>180 & bestFitAngs2D_AllData<=270)-180;
bestFitAngs2D_AllData(bestFitAngs2D_AllData>270 & bestFitAngs2D_AllData<=360)...
    =bestFitAngs2D_AllData(bestFitAngs2D_AllData>270 & bestFitAngs2D_AllData<=360)-270;



% Plot histogram
figure(fig); clf; hold on; fig=fig+1;
plotTitle = 'Population 2-D Fit to All Data';
set(gcf,'Name',plotTitle,'NumberTitle','off')
hist(bestFitAngs2D_AllData,32);
xlabel('Best Fit Rotation (deg)')
ylabel('Number of Cells')


%% Cone Weight Analysis

% for d = 1:size(data,1)
%     
%     datafile = data{d}{1}.datafile;
%     disp('ConeWeightAnalysis')
%     disp(['Now processing ' datafile])
%         
%     data{d}{3} = ConeWeightAnalysis(data{d}{3});
%     
%     for p = 1:numel(data{d}{3})
%         if d==1 && p==1
%             CWAPos = data{d}{3}{p}.ConeWeightAnal.Pos(1,:);
%             CWANeg = data{d}{3}{p}.ConeWeightAnal.Neg(1,:);
%         end
%         CWAPos = cat(1,CWAPos,data{d}{3}{p}.ConeWeightAnal.Pos(end,:));
%         CWANeg = cat(1,CWANeg,data{d}{3}{p}.ConeWeightAnal.Neg(end,:));
%     end
% 
% end
% 
% % Some variables for plotting
% maxPos = max(max(cell2mat(CWAPos(2:end,3))),max(cell2mat(CWAPos(2:end,4))));
% maxNeg = max(max(cell2mat(CWANeg(2:end,4))),max(cell2mat(CWANeg(2:end,4))));
% 
% 
% % Plot Results
% figure(fig); clf; hold on; grid on; fig=fig+1;
% plotTitle = 'Population Cone Weight Analysis';
% set(gcf,'Name',plotTitle,'NumberTitle','off')
% subplot(1,2,1); hold on; grid on;
% plot(cell2mat(CWAPos(2:end,3)),cell2mat(CWAPos(2:end,4)),'*')
% plot(0:.1:maxPos,0:.1:maxPos,'k--')
% xlabel('Response to L-Cone Isolating Stimuli (FR)')
% ylabel('Response to M-Cone Isolating Stimuli (FR)')
% title('Positive Cone-Isoating Contrasts')
% subplot(1,2,2); hold on; grid on;
% plot(cell2mat(CWANeg(2:end,3)),cell2mat(CWANeg(2:end,4)),'*')
% plot(0:.1:maxNeg,0:.1:maxNeg,'k--')
% xlabel('Response to L-Cone Isolating Stimuli (FR)')
% ylabel('Response to M-Cone Isolating Stimuli (FR)')
% title('Negative Cone-Isoating Contrasts')
 

%% Second Cone Weight Analysis
% 
% mostSensitiveDir_AllData = [];
% 
% for d = 1:size(data,1)
%     
%     datafile = data{d}{1}.datafile;
%     disp('2-D Fit To All Data')
%     disp(['Now processing ' datafile])
%                 
%     data{d}{3} = ConeWeightAnalysis2(data{d}{3});
%     
%     for p = 1:numel(data{d}{3})
%         
%         tempAng = data{d}{3}{p}.ConeWeightAnal2.MostSensitiveDirInDeg;
%         
%         mostSensitiveDir_AllData = cat(1,mostSensitiveDir_AllData,tempAng);
%     
%     end
% 
% end
% 
% figure(fig); clf; hold on; grid on; fig=fig+1;
% plotTitle = 'Population Most Sensitive Direction';
% set(gcf,'Name',plotTitle,'NumberTitle','off')
% %hist(mostSensitiveDir_AllData,floor(size(data,1)/2))
% hist(mostSensitiveDir_AllData,0:22.5:360)
% xlabel('Most Sensitive Direction (Deg)')
% ylabel('Number of Cells')


%% 2-D Fit to Two Orthogonal Axes
% 
% bestFitAngs_TwoAxes = [];
% 
% for d = 1:size(data,1)
%     
%     datafile = data{d}{1}.datafile;
%     disp('2-D Fit To Two Orthogonal Axes')
%     disp(['Now processing ' datafile])
%                 
%     data{d}{3} = TwoDFit_TwoAxes(data{d}{3});
%     
%     [LL maxidx] = max(cell2mat(data{d}{3}{1}.TwoDFitTwoAxes(2:end,2)));
%     tempAng = cell2mat(data{d}{3}{1}.TwoDFitTwoAxes(maxidx+1,1));
%     
%     bestFitAngs_TwoAxes = cat(1,bestFitAngs_TwoAxes,tempAng);
% 
% end
% 
% figure(fig); clf; hold on; grid on; fig=fig+1;
% plotTitle = 'Population 2-D Fit to Two Orthogonal Axes';
% set(gcf,'Name',plotTitle,'NumberTitle','off')
% hist(bestFitAngs_TwoAxes)
% xlabel('Best Fit Rotation (deg)')
% ylabel('Number of Cells')
%     