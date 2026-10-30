function Main_particleAnalysis(particleExample_or_cellExample,WindowsOrMac)

%% Section 1: Get Initial File List, configure paths, configure initial directoires.

currentFolder = pwd;

% location of example data
% windows=1, mac=0, this to correct for '/'
if WindowsOrMac==0
    inputDir = [currentFolder '/data'];
    outputDir = [currentFolder '/results'];
    addpath(genpath([currentFolder '/data']))
elseif WindowsOrMac==1
    inputDir = [currentFolder '\data'];
    outputDir = [currentFolder '\results'];
    addpath(genpath([currentFolder '\data']))
end

% location of where to the processed data will be saved

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% find the data in the folder, here the data is saved as .mat files
% particleExample_or_cellExample=1;
if particleExample_or_cellExample==1
    fileType = 'partilceData*';  % particles on glassin water
elseif particleExample_or_cellExample==2
    fileType = 'cellParticleData*'; % particles in cells
end

fileList = dir(fullfile(inputDir, fileType));

% Section 2: This Section tracks files in data directory, saves mat into outputDir
for k =  1% 1 :length(fileList)
    % Get the full path of the current .nd2 file
    inputFile = fullfile(inputDir, fileList(k).name);
    % Extract the base file name without extension
    [~, baseFileName, ~] = fileparts(fileList(k).name);
    outputFileName = [baseFileName '_tracked'];
    automateTrack_particle_Fluo_DIC_iSCATandDF_clean(inputFile, outputFileName, outputDir,WindowsOrMac);

end


%% once the data is processed, we can create different plots
WindowsOrMac=0; % windows=1, mac=0, this to correct for '/'
% particleExample_or_cellExample=1;

% thresholds for the goodness of fit for the included particles
Field_GOF_thres=0.6;



if particleExample_or_cellExample==1
    if WindowsOrMac==1
        load([currentFolder '\results\partilceData_100 nm polystryene on glass_tracked.mat'])
    else
        load([currentFolder '/results/partilceData_100 nm polystryene on glass_tracked.mat'])
    end
    DIC_GOF_thres=0.4; % if particle on glass analysis

elseif particleExample_or_cellExample==2
    if WindowsOrMac==1
        load([currentFolder '\results\cellParticleData_100nmPolystyrene_tracked.mat'])
    else
        load([currentFolder '/results/cellParticleData_100nmPolystyrene_tracked.mat'])
    end
    DIC_GOF_thres=-1; % if particle in cell analysis
end

% code that sets the z=0 position based on the maximum background signal
for jj=1:size(bg_vec2,1)
    test_x=1:length(bg_vec2(jj,:));
    spline_x=linspace(1,max(test_x),10*length(test_x)-1);
    s = spline(test_x,bg_vec2(jj,:),spline_x);
    [a,b]=max(s);
%     sel=Positions(:,10)==jj;
    Positions(:,8)=(Positions(:,8)-spline_x(b))*px_size_Z;
end

% the array with the data is formatted as 
% index (1,2,8) corresponds to (x,y,z)
% index (12,13) corresponds to C-SCAT estimates (405 and 660 nm)
% index (3,4) corresponds to iSCAT estimates (405 and 660 nm)
% index (5,6) corresponds to DIC and fluoresence
% index (7,14) are DIC and C-SCAT goodness of fit estiamtes
% index (9-11) local background, video number (in case of multiple videos), signal to background ratio



sel_low=((Positions(:,13)>1).*(Positions(:,8)>=1.0).*(Positions(:,14)>=Field_GOF_thres).*(Positions(:,6)>10).*(Positions(:,7)>DIC_GOF_thres))==1;
sel_high=((Positions(:,13)>1).*(Positions(:,8)<1.0).*(Positions(:,14)>=Field_GOF_thres).*(Positions(:,6)>10).*(Positions(:,7)>DIC_GOF_thres))==1;
sel=(1*sel_low+sel_high)>0;
both_populations_present=(sum(sel_low).*sum(sel_high))>0;


% plot C-SCAT vs z-position, useful for cell data
figure(1); clf; hold on
if both_populations_present==0
    scatter(2*Positions(sel,13)/sqrt(1000),Positions(sel,8),20,'r','o','filled', 'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.0)
else
    scatter(2*Positions(sel_low,13)/sqrt(1000),Positions(sel_low,8),20,'b','o','filled', 'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.0)
    scatter(2*Positions(sel_high,13)/sqrt(1000),Positions(sel_high,8),20,'r','o','filled', 'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.0)
    legend('z>1 \mum','z\leq1 \mum')
end
ylabel('Z-position [\mum]','Interpreter','tex')
xlabel('Integrated C-SCAT signal','Interpreter','tex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
box on
set(gca,'linewidth',2)
set(gca,'FontSize', 15)
set(gcf, 'PaperPositionMode', 'auto');
grid off

xlim([0,80])
ylim([-0.5,8])

% if particleExample_or_cellExample==1
%     saveas(gcf, '../results/particleOnGlass_Cscat_vs_z.png');
% else
%     saveas(gcf, '../results/particleInCells_Cscat_vs_z.png');
% end

%
% plot the signal to background ratio
hist_vec=logspace(-1.8,2,20);
figure(2); clf; hold on

if both_populations_present==0
    h1=histogram(Positions(sel,11),hist_vec,'Normalization','probability');
    h1.EdgeColor = 'none';

else

    h1=histogram(Positions(sel_low,11),hist_vec,'Normalization','probability');
    h1.EdgeColor = 'none';

    h2=histogram(Positions(sel_high,11),hist_vec,'Normalization','probability');
    h2.EdgeColor = 'none';
    legend('z>1 \mum','z\leq1 \mum')
end


set(gca,'Xscale','log') 
ylabel('Probability','Interpreter','tex')
xlabel('|s|/r','Interpreter','tex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(gca,'linewidth',2)
set(gca,'FontSize', 15)
set(gcf, 'PaperPositionMode', 'auto');
grid off
box on


xlim([min(hist_vec),max(hist_vec)])

% if particleExample_or_cellExample==1
%     saveas(gcf, '../results/particleOnGlass_signal_to_bg_ratio.png');
% else
%     saveas(gcf, '../results/particleInCells_signal_to_bg_ratio.png');
% end

% plot the C-SCAT signal histogram
hist_vec=0:3:100;
figure(3); clf; hold on

if both_populations_present==0
    h2=histogram(2*Positions(sel,13)/sqrt(1000),hist_vec,'Normalization','probability');
    h2.EdgeColor = 'none';
else
    h1=histogram(2*Positions(sel_low,13)/sqrt(1000),hist_vec,'Normalization','probability');
    h1.EdgeColor = 'none';

    h2=histogram(2*Positions(sel_high,13)/sqrt(1000),hist_vec,'Normalization','probability');
    h2.EdgeColor = 'none';

    legend('z>1 \mum','z\leq1 \mum')
end


ylabel('Normalized occurance','Interpreter','tex')
xlabel('Integrated particle signal','Interpreter','tex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
box on
set(gca,'linewidth',2)
set(gca,'FontSize', 15)
set(gcf, 'PaperPositionMode', 'auto');
grid off


ylim([0,0.3])
xlim([0,80])


% if particleExample_or_cellExample==1
%     saveas(gcf, '../results/particleOnGlass_Cscat_histogram.png');
% else
%     saveas(gcf, '../results/particleInCells_Cscat_histogram.png');
% end
%

% plot C-SCAT vs background amplitude
figure(4); clf; hold on

if both_populations_present==0
    scatter(2*Positions(sel,13)/sqrt(1000),Positions(sel,9),50,'o', 'MarkerFaceColor','r', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
else
    scatter(2*Positions(sel_low,13)/sqrt(1000),Positions(sel_low,9),50,'o', 'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
    scatter(2*Positions(sel_high,13)/sqrt(1000),Positions(sel_high,9),50,'o', 'MarkerFaceColor','r', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
    legend('z>1 \mum','z\leq1 \mum')
end


xlim([0,80])
ylim([1,30000])

ylabel('Background amplitude','Interpreter','tex')
xlabel('Integrated C-SCAT signal','Interpreter','tex')
set(gca, 'YScale', 'log')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
box on
set(gca,'linewidth',2)
set(gca,'FontSize', 15)
set(gcf, 'PaperPositionMode', 'auto');
grid off
box on


% if particleExample_or_cellExample==1
%     saveas(gcf, '../results/particleOnGlass_Cscat_vs_bg.png');
% else
%     saveas(gcf, '../results/particleInCells_Cscat_vs_bg.png');
% end

% plot C-SCAT vs DIC
figure(5); clf; hold on

if both_populations_present==0
    scatter(2*Positions(sel,13)/sqrt(1000),Positions(sel,6),20,'o','filled', 'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.0)
else
    scatter(2*Positions(sel_low,13)/sqrt(1000),Positions(sel_low,6),50,'o', 'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
    scatter(2*Positions(sel_high,13)/sqrt(1000),Positions(sel_high,6),50,'o', 'MarkerFaceColor','r', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.0)
    legend('z>1 \mum','z\leq1 \mum')
end


xlim([0,80])
ylim([0,400000])

ylabel('Fluoresence','Interpreter','tex')
xlabel('Integrated C-SCAT signal','Interpreter','tex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
box on
set(gca,'linewidth',2)
set(gca,'FontSize', 15)
set(gcf, 'PaperPositionMode', 'auto');
grid off

% if particleExample_or_cellExample==1
%     saveas(gcf, '../results/particleOnGlass_Cscat_vs_fluo.png');
% else
%     saveas(gcf, '../results/particleInCells_Cscat_vs_fluo.png');
% end

end