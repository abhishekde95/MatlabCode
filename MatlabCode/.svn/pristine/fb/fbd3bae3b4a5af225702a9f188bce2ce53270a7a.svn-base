% ForShane
% Based on FACSworkup.m
% GDLH 12/20/19

folder_name = '/Users/horwitzlab/Desktop/NAb Testing';
old_dir = pwd;
disp(['Currently in ',old_dir])
cd(folder_name);
[fname_nc,pathname] = uigetfile('*.*','Which file for negative control?','MultiSelect','off');
if fname_nc == 0
   error('No negative control file selected. Aborting.'); 
end
Data = importdata([pathname,filesep,fname_nc]);
cd(old_dir);

% Gating the events
FSCA = Data.data(:,strcmp(Data.textdata,'"FSC-A"'));
FSCH = Data.data(:,strcmp(Data.textdata,'"FSC-H"'));
FSCW = Data.data(:,strcmp(Data.textdata,'"FSC-W"'));
SSCA = Data.data(:,strcmp(Data.textdata,'"SSC-A"'));
GFP = Data.data(:,strcmp(Data.textdata,'"GFP-A"'));
DAPI = Data.data(:,strcmp(Data.textdata,'"DAPI-A"'));
MCHERRY = Data.data(:,strcmp(Data.textdata,'"mCherry-A"'));
TIME = Data.data(:,strcmp(Data.textdata,'"TIME"'));


%% Gate 1
h_fig = figure;
hold on; box off
plot(FSCA,SSCA,'.k','MarkerSize',1)
set(gca,'TickDir','out')
title('Draw a Polygon around the Cluster of Interest')
xlabel('Forward Scatter')
ylabel('Side Scatter');
btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
        'Position', [20 20 60 20],...
        'Callback', 'uiresume');
h = impoly; %Allows the user to interactively draw a polygon on the data.

%Pausing the analysis until the gate is completed and the user decides to
%continue
uiwait

%Calling the vertices of the polygon 'Gate1'
Gate1 = getPosition(h);

%Defining the x and y components of the Gate
xv = Gate1(:,1);
yv = Gate1(:,2);

%inpolygon takes inputs of the x and y data as well as the x and y
%positions of the gate to output a logical vector which indicates whether a
%particular point was included or excluded from the center of the polygon
%defined by the gate.
Lgate1 = inpolygon(FSCA,SSCA,xv,yv);

%% Gate 2
%Plotting a new subspace of the reduced data set
figure(h_fig);
clf;
MyClr = [0,0.8,0.6];
hold on; box off
plot(FSCH(Lgate1),FSCW(Lgate1),'.k','MarkerSize',1.5)
set(gca,'TickDir','out')
title('Draw a Polygon around the Cluster of Interest')
xlabel('Forward Scatter Height')
ylabel('Forward Scatter Width')

%Creating a bottum that will resume the analysis after being pressed
btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
        'Position', [20 20 60 20],...
        'Callback', 'uiresume');
    
%Drawing a new polygon
i = imrect;

%Pausing the analysis until the gate is completed and the user decides to
%continue
uiwait
%Establishing the coordinates for gate 2 given from impoly
Gate2 = getPosition(i);
Gate2 = [Gate2(1) Gate2(2)
         Gate2(1)+Gate2(3) Gate2(2)
         Gate2(1)+Gate2(3) Gate2(2)+Gate2(4)
         Gate2(1) Gate2(2)+Gate2(4)];
     
% Lgate_both = inpolygon(FSCH(Lgate1),FSCW(Lgate1),Gate2(:,1),Gate2(:,2));
Lgate_both = Lgate1 & inpolygon(FSCH,FSCW,Gate2(:,1),Gate2(:,2));


%%
% Which reporter signal
reporteridx = listdlg('PromptString','Which reporter?',...
                      'SelectionMode','single',...
                      'ListString',{'GFP','mCherry'});
if reporteridx == 1
    reporter_data = GFP;
    reporter_label = 'GFP';
else
    reporter_data = MCHERRY;
    reporter_label = 'mCherry';
end
%% 
% Plotting the reporter signal for the negative control
figure(h_fig);
clf;

%L_sensible_value = reporter_data > 0;
%bins = linspace(min(log10(reporter_data(L_sensible_value & Lgate_both)))-.5, max(log10(reporter_data(L_sensible_value & Lgate_both)))+.5,100);
%[n1,x] = hist(log10(reporter_data(L_sensible_value & Lgate_both)),bins);
if any( reporter_data < 0 )
    error('cannot have negative values');
end

bins = linspace(min(log10(reporter_data(Lgate_both)))-.5, max(log10(reporter_data(Lgate_both)))+.5,100);
[n1,x] = hist(log10(reporter_data(Lgate_both)),bins);

bar(x,n1);
set(gca,'Xlim',[bins(1) bins(end)],'Yscale','linear');
[threshold, ~] = ginput(1);
neg_cntrl_counts_and_total = [sum(log10(reporter_data(Lgate_both)) > threshold) sum(Lgate_both)];
xlabel([reporter_label,' signal']);
ylabel('counts');
title(['Negative control: ',fname_nc]);
%%
% Loading experimental data
cd(folder_name);
[fnames,pathname] = uigetfile('*.*','Which experimental files?','MultiSelect','on');
cd(old_dir);

condition_names = {};
for i = 1:length(fnames)
   condition_names(i) = inputdlg(['Enter condition for ',fnames{i}]);
end

figure(h_fig);
set(gcf,'Position',[128    90   440   890])
data = neg_cntrl_counts_and_total;
for i = 1:length(fnames)
    fname  = fnames{i};
    Data = importdata([pathname,filesep,fname]);
    FSCA = Data.data(:,strcmp(Data.textdata,'"FSC-A"'));
    FSCH = Data.data(:,strcmp(Data.textdata,'"FSC-H"'));
    FSCW = Data.data(:,strcmp(Data.textdata,'"FSC-W"'));
    SSCA = Data.data(:,strcmp(Data.textdata,'"SSC-A"'));
    GFP = Data.data(:,strcmp(Data.textdata,'"GFP-A"'));
    MCHERRY = Data.data(:,strcmp(Data.textdata,'"MCHERRY-A"'));

    Lgate1 = inpolygon(FSCA,SSCA,xv,yv);
    Lgate_both = Lgate1 & inpolygon(FSCH,FSCW,Gate2(:,1),Gate2(:,2));
    subplot(length(fnames),1,i); hold on;
    if reporteridx == 1
        reporter_data = GFP;
    else
        reporter_data = MCHERRY;        
    end
    if ~isempty(reporter_data)
        [n2,x] = hist(log10(reporter_data(Lgate_both)),bins);
        bar(x,n2);
        plot(threshold, max(n2),'kv');
        set(gca,'Xlim',[bins(1) bins(end)]);
        data = [data; sum(log10(reporter_data(Lgate_both)) > threshold) sum(Lgate_both)];
        title(condition_names{i});
        ylabel('counts');
        set(gca,'Yscale','log');
    end
end
equatesubplotaxeslims;
xlabel([reporter_label,' signal']);
figure; axes; hold on;
plot(data(:,1)./data(:,2),'o','MarkerFace','black');
ylabel('Proportion GFP+');
xlabel('Condition');
set(gca,'Xtick',1:length(fnames)+1,'XtickLabel',cat(2,'NC',condition_names));


