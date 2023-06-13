%% Load environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/6/21 DAC                                    %   
% 9/6/21 modified by MDF                        %
% 10/11/22 modified by JEK from Unfold2021_v3.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
who_is_running = 'J'; % 'M': Maria or Damian or Juan % jek 10/11/22
sanity_checks = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(who_is_running,'M')
    cd  /home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/functions
elseif strcmp(who_is_running,'D')
    cd /Users/dac/Documents/GitLab/corregistro/Analysis_dac/functions   % mdf 8/6/21 DAC directory
elseif strcmp(who_is_running,'J')
    cd /home/juank/Desktop/EJN2022/matlab/                              % jek 10/11/22
end

[runpath,code_path,session_path,whoisrunning]=add_paths_EEG_ET_beta(0,'J');% edited by jek 10/11/22
addpath(code_path.eeglabpath); 
eeglab
init_unfold

addpath(code_path.toolbox); 
addpath(genpath(code_path.corregistro_fn));
addpath(code_path.my_functions);
addpath(code_path.psyco);
addpath(runpath);
addpath(code_path.functions);

clear fp; fp = FP_epochAnalysis();

matdir              = session_path.matfiles;
datapath            = session_path.data_analysis; 
% auxdir              =  session_path.aux_data; % mdf 8/6/21 DAC directory

if strcmp(who_is_running,'M');      addpath(genpath('/home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/UNFOLD'))
elseif strcmp(who_is_running,'D');  addpath(genpath('/Users/dac/Documents/GitLab/corregistro/Analysis_dac/UNFOLD')) % mdf 8/6/21 DAC directory
elseif strcmp(who_is_running,'J');  addpath(genpath('/home/juank/repos/corregistro/Analysis_dac/UNFOLD'))           % jek 10/11/22
end

fp.cfg.subjname         = session_path.subjname;
fp.cfg.sessionfilenames = session_path.sessionfilenames;
fp.cfg.matdir           = matdir;
fprintf('****************************************************************\n');    
%Edit this lines to choose a model to fit It should be written in Wilkinson
%notation, categorical and variables model with splines should be indicated
%cat(categorical variable), spl(variable,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fig. 6
% % fixs_csv_file  = '16subj_VSNTpostvsEX_for_unfold_GlobalStandarisedCenteredRank.csv';% OLD This is the file used in EJN Round 0 for VSNTpre_EX condition. It has 5127 fixs.
% % fixs_csv_file  = 'EJN2022_16subj_VSNT_EX.csv';                                      % NEW This is the file used in EJN Round 1 for VSNT_EX condition. It has 6421 fixs.
% % fixs_csv_file  = 'EJN2022_16subj_VSNTpre_EX.csv';                                   % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% fixs_csv_file  = 'EJN2022_16subj_VSNTpre_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'pre_';
%
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+faces+NCrank+interTaskNCRank';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank
% 
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% variablesOI     = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'};
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% cond            = {'EX','AS'};%%%cond = split(hipothesis,'_');

% % New analysis: Interaction between Task and Category (same data)
% fixs_csv_file   = 'EJN2022_16subj_VSNTpre_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'pre_';
% 
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+faces+interTaskFaces';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + faces:VSNT_EX'};
% 
% variablesOI     = {'(Intercept)' hipothesis 'faces' 'VSNT_EX:faces'};
% cond            = {'EX','AS'};%%%cond = split(hipothesis,'_');

% New analysis: Full model like figure 6 but with absent trials instead pretarget
fixs_csv_file   = 'EJN2022_16subj_VSNTpre_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
datasetID       = 'abs_';

hipothesis      = 'VSNT_EX';
modelName       = 'Inter+task+faces+NCrank+interTaskNCRank';
model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank

%e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
variablesOI     = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'};
%e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
cond            = {'EX','AS'};%%%cond = split(hipothesis,'_');

fprintf('CSV file: %s\n', fixs_csv_file);
fprintf('Models name: %s\n', modelName);
fprintf('Model (Wilkinson notation): %s\n', model{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baselineCorr    = 1;
% ourPermTest     = 1;
fprintf('Analysing with baseline correction\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Choose the freq-band
Theta           = struct('limits',[4  8],'name','Theta');
Alpha           = struct('limits',[8 13],'name','Alpha');
Beta            = struct('limits',[15 25],'name','Beta');
LowBeta         = struct('limits',[14 19],'name','LowBeta');
HighBeta        = struct('limits',[20 30],'name','HighBeta');
LowGamma        = struct('limits',[30 40],'name','LowGamma');

currentBand     = Theta;

freqFilteredEEG = 0; % 1 if we used the EEG filtered in freq

if freqFilteredEEG
    fprintf('Analysing %s band\n',currentBand)
    inFolder    = [session_path.freqEEG currentBand.name '/'];
    names       = fp.loadSubjects(inFolder, 1, 1, 'set');
    subj_index  = [1:length(names)];
else   
    fprintf('Analysing broad band (FRPs)\n')
    inFolder    =  session_path.data_analysis;
    names       = fp.loadSubjects(inFolder, 1, 1, 'set');
    subj_index  = [1:4 6:11 13:15 19 22:23];
end
Nsuj = length(subj_index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter conditions % jek 10/11/2022
filter_faces    = 0; % 1 if filter fix to faces O otherwise
filter_EX       = 0;
filter_AS       = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jek 10/11/2022: Auxiliary files??? These are not necessarily the same
if strcmp(who_is_running,'M');      auxdir      =  session_path.fixs;           % fixs
elseif strcmp(who_is_running,'D');  auxdir      =  session_path.aux_data;       % Analysis_aux_data
elseif strcmp(who_is_running,'J');  auxdir      =  session_path.fixs;           % python: the csv with the output of select_fixs_unfold-EJN2022.ipynb
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main output folder
folderUnf               = fullfile(session_path.out, 'Unfold2022/' );           % jek 10/11/22: Partial output. I created this folder.
    if ~exist(folderUnf,'dir');     mkdir(folderUnf); end

% EEG data folder for UNFOLD (1st STEP, it needs to run once per dataset)    
if freqFilteredEEG; folderPrepared      = fullfile(folderUnf, [datasetID, 'preparedEEG/'], currentBand.name, '/');
else;               folderPrepared      = fullfile(folderUnf, [datasetID, 'preparedEEG/'] );
end
    if ~exist(folderPrepared,'dir');mkdir(folderPrepared); end

% Models output folder for UNFOLD (remaining STEPs)        
models = fullfile(folderUnf, [datasetID, 'models/'] );
    if ~exist(models, 'dir'); mkdir(models); end

if freqFilteredEEG; folderModelName     = [currentBand.name '/' modelName '/'];
else;               folderModelName     = [modelName '/'];
end

folderOut               = fullfile(models,folderModelName);
    if ~exist(folderOut, 'dir'); mkdir(folderOut); end

fprintf('****************************************************************\n');

%% Heatmaps
fileStruct      = fullfile(folderOut,'variableStruct.mat');

load(fileStruct);
times   = variableStruct.times; 
%     chanlocs = load(fullfile(session_path.aux_data,'locationFile.mat'));
chanlocs = load('locationFile.mat');
chanlocs = chanlocs.e_loc;

variables       = fieldnames(variableStruct);
variables(cellfun(@(x) strcmp(x,'times'),variables)) = [];
variables(cellfun(@(x) strcmp(x,'Suj'),variables)) = [];

% limits of significant time-window defined by visual inspection
xmin_sig = 0.0;
xmax_sig = 0.0;
    
for ii = 1:length(variables)
    variable2TFCE   = variables{ii};

    TFCEfilename = ['tfceResults_' variable2TFCE '.mat' ];
    fileTFCE  = fullfile(folderOut, TFCEfilename);
    load(fileTFCE)


    Y = tfceResult.P_Values;
    figure(10 + ii); 
        imagesc(times,1:64,Y([1:32 64:-1:33],:),[0 0.05]); 
        
        % Colormap
        colormap hot; 
        c=colorbar; 
        ylabel(c,'p value','fontsize',42); 
        set(c,'fontsize',35)
    
        xlabel('Time [s]','fontsize',42)
        ylabel('Channels','fontsize',42)
        title(variable2TFCE,'fontsize',42)
        yticks([7 20 27 33 40 55])
        nuChans = {chanlocs([1:32 64:-1:33]).labels};
        yticklabels(nuChans([7 25 27 33 40 62]))
        set(gcf,'Color','w')
        set(gca,'fontsize',35)
        xlim([-0.1,0.4])
        xline(xmin_sig,'Color',[0 0 0.5]+0.005,'linewidth', 5)
        xline(0,'-k','linewidth', 5)
        xline(xmax_sig,'Color',[0 0 0.5]+0.005,'linewidth', 5)
        ax = gca;
        get(gca)
        set(gca,'linewidth', 5)

end


%% topoplot significant channels

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_window = [xmin_sig xmax_sig];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inFolder                =  fullfile(session_path.out, 'Unfold2022/' );
models = fullfile(inFolder, 'models/' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use the name of the folder where the model were stored
model2analyse = modelName;

if freqFilteredEEG
    folderOut = fullfile(models,currentBand.name, '/',model2analyse);
else   
    folderOut = fullfile(models,model2analyse);
end

fileStruct  = fullfile(folderOut,'variableStruct.mat');
if exist(folderOut, 'dir')
    load(fileStruct);
    times   = variableStruct.times; 
%     chanlocs = load(fullfile(session_path.aux_data,'locationFile.mat'));
    chanlocs = load('locationFile.mat');
    chanlocs = chanlocs.e_loc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if freqFilteredEEG
    folderIn = fullfile(models,currentBand.name, '/',model2analyse);
else   
    folderIn = fullfile(models,model2analyse);
end
TFCEfilename = ['tfceResults_' variable2TFCE '.mat' ];
fileTFCE  = fullfile(folderIn, TFCEfilename);
load(fileTFCE)

y=sum(tfceResult.P_Values(:,times>time_window(1) & times<time_window(2))<0.05,2);

y = y/max(y);

figure; topoplot(y,chanlocs,'maplimits',[0 1]); colormap('gray'); %colorbar;

%tiledlayout(1,1,'Padding','tight');

%% plot invividual FRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use the name of the folder where the model were stored
model2analyse = modelName;
inFolder                =  fullfile(session_path.out, 'Unfold2022/' );
models = fullfile(inFolder, 'models/' );
if freqFilteredEEG
    folderIn = fullfile(models,currentBand.name, '/',model2analyse);
else   
    folderIn = fullfile(models,model2analyse);
end
fileStruct  = fullfile(folderIn,'variableStruct.mat');

if exist(folderIn, 'dir')
    load(fileStruct);
    times   = variableStruct.times; 
    chanlocs = load('locationFile.mat');
%     chanlocs = load(fullfile(session_path.aux_data,'locationFile.mat'));
    chanlocs = chanlocs.e_loc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variable2plot =  {'Intercept' 'faces' hipothesis 'NCRank' 'NCRank_VSNT_EX'} %'StandFixRank_VSNT_EX' 'StandFixRank_isEX'}

TFCEbar = 0; % to plot the TFCE result of the last variable included
if TFCEbar
    if sum(strcmp(variable2plot,'Intercept'))>0
        tfce = cell(length(variable2plot)-1,1);
    else
        tfce = cell(length(variable2plot),1);
    end
    for var = 1:length(variable2plot)
        if ~strcmp(variable2plot{var},'Intercept')
            TFCEfilename = ['tfceResults_' variable2plot{var} '.mat' ];
            fileTFCE  = fullfile(folderIn, TFCEfilename);
            load(fileTFCE)
            if sum(strcmp(variable2plot,'Intercept'))>0
                tfce{var-1} = tfceResult;
            else
                tfce{var} = tfceResult;
            end
        end
    end
end

chan = {'PO7','Fz','PO8','TP7','CPz','TP8','PO7','Pz','PO8','O1','Oz','O2','F1','F2','Cz','AFz','POz'};

f=fields(variableStruct);

XLIMI = [-.100 .400];
% yloglim = [10^-4  10^0];
set(gcf,'Color','w')
erppos = {[ 1 3 5],[2 4 6],[9 11 13],[10 12 14]} ;

if strcmp(hipothesis,'faces') || strcmp(hipothesis,'NTF_NTO')
    chanplot = [15 2 1 9 10 12] %NTFvsNTO
elseif strcmp(hipothesis,'VSNT_EX') || strcmp(hipothesis,'VSNTpre_VSNTpost') || strcmp(hipothesis,'VSNTpost_EX') || strcmp(hipothesis,'isEX')
    chanplot = [15 2 13 14 10 12 16 17 11] % {chanlocs.labels} % ASvsEX
elseif strcmp(hipothesis,'target_AS')
    chanplot = [2 15 8 11];
end
 %13 14 10 12  ]

for ifig = 1:9
    figN = 102 + ifig;
    figure(figN); clf
    i = chanplot(ifig)
    ch = find(ismember({chanlocs.labels},chan{i}));
    
    for var = 1:length(variable2plot)
        struct(var).x = []
    end
    for iSuj = [1:16]
        for var = 1:length(variable2plot)
            struct(var).x = [struct(var).x; variableStruct(iSuj).(variable2plot{var})(ch,:)];
        end
    end
    
    for var = 1:length(variable2plot)
        struct(var).X = squeeze(mean(struct(var).x,1));
        struct(var).eX = squeeze(std(struct(var).x,0,1))/sqrt(size(struct(var).x,1));
        struct(var).base = mean(struct(var).X(-.2<times & times<0));
        if baselineCorr==1
            struct(var).X = struct(var).X - struct(var).base
        end
        chan{i}
    end   

    ylabel(strcat(chan{i},' [\muV]'),'FontSize',22) 

    box on                                    
    hold on
    %title(chan{i})   
    
    for var = 1:length(variable2plot)
        struct(var).ylim = max([abs(max(abs(struct(var).X-struct(var).eX),[],'all')),abs(max(abs(struct(var).X+struct(var).eX),[],'all'))]);
    end
 
    ylim_ = max([struct.ylim])
    YLIMI = [-ylim_-ylim_/2 ylim_];
    
    colors = {'k-' 'r-' 'g-' 'b-' 'y-'};
    x = [times,     times(end:-1:1),            times(1)];
    for var = 1:length(variable2plot)
        y = [struct(var).X+struct(var).eX,    struct(var).X(end:-1:1)-struct(var).eX(end:-1:1), struct(var).X(1)+struct(var).eX(1)];
        patch(x,y,colors{var},'EdgeColor','none','FaceAlpha',0.5)
    end

    h = [];
    hlegend = [];
    for var = 1:length(variable2plot)
        h(var) = plot(times,struct(var).X,colors{var},'LineWidth',3,'DisplayName',variable2plot{var});
        hlegend = [hlegend h(var)];
        if TFCEbar
            if ~strcmp(variable2plot{var},'Intercept')
                if sum(strcmp(variable2plot,'Intercept'))>0
                    p_block = find([tfce{var-1,1}.P_Values(ch,:)]<0.05);
                else
                    p_block = find([tfce{var,1}.P_Values(ch,:)]<0.05);
                end
                difp = diff(p_block);
                indcut=1;
                for j = 1:length(difp)
                    if difp(j) ~= 1
                        indcut = [indcut j];
                    end
                end
                indcut = [indcut length(difp)];
                if sum(strcmp(variable2plot,'Intercept'))>0
                    Y = tfce{var-1,1}.P_Values(ch,:); Y(Y>0.05)=1;
                else
                    Y = tfce{var,1}.P_Values(ch,:); Y(Y>0.05)=1;
                end
                imagesc(times(-99<times & times<1000*XLIMI(2)-1),[-5.5 -5.3],Y(find(-99<times & times<1000*XLIMI(2)-1)),[0 0.05]); colormap hot;
                if ~freqFilteredEEG
                    imagesc(times(-99<times & times<1000*XLIMI(2)-1),[-5 -4.5],Y(find(-99<times & times<1000*XLIMI(2)-1)),[0 0.05]); colormap hot;
                end
            end
        end
    end

    ax=gca;
    ax.FontSize = 18; 

    hold off

    a = get(gca,'XTickLabel');
    xlabel('time [s]','FontSize',18);
    xlim([-0.1,0.4])
    ylim([-6,6])
    xline(0,'-k','linewidth', 3)
    yline(0,'-k','linewidth', 3)
    set(gca, 'PlotBoxAspectRatio', [1,1,1],'linewidth', 3)
    get(gca)

%     saveas(gca,strcat('/media/cbclab/MARIADAFON1T/Analysis_2020/Unfold2021/models/',currentBand.name,'/Inter+task+faces+NCrank+interTaskNCRank/',chan{i},'_FRP.jpg'));
    saveas(gca,strcat('/home/juank/Desktop/EEGEYE_EJN2022/Unfold2022/models/Inter+task+faces+NCrank+interTaskNCRank/',chan{i},'_FRP.jpg'));
end    

%% plot an arrange of FRPs (old-code)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use the name of the folder where the model were stored
model2analyse = modelName;
inFolder                =  fullfile(session_path.out, 'Unfold2021/' );
models = fullfile(inFolder, 'models/' );
if freqFilteredEEG
    folderIn = fullfile(models,currentBand.name, '/',model2analyse);
else   
    folderIn = fullfile(models,model2analyse);
end
fileStruct  = fullfile(folderIn,'variableStruct.mat');

if exist(folderIn, 'dir')
    load(fileStruct);
    times   = variableStruct.times; 
    chanlocs = load(fullfile(session_path.aux_data,'locationFile.mat'));
    chanlocs = chanlocs.e_loc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variable2plot = {'Intercept' 'faces'};

if sum(strcmp(variable2plot,'Intercept'))>0
    tfce = cell(length(variable2plot)-1,1);
else
    tfce = cell(length(variable2plot),1);
end
for var = 1:length(variable2plot)
    if ~strcmp(variable2plot{var},'Intercept')
        TFCEfilename = ['tfceResults_' variable2plot{var} '.mat' ];
        fileTFCE  = fullfile(folderIn, TFCEfilename);
        load(fileTFCE)
        if sum(strcmp(variable2plot,'Intercept'))>0
            tfce{var-1} = tfceResult;
        else
            tfce{var} = tfceResult;
        end
    end
end

ventanas={[0.07 0.12],[0.12 0.15],[0.15 0.2],[0.2 0.23]};      

figure(102); clf

chan = {'PO7','Fz','PO8','TP7','CPz','TP8','PO7','Pz','PO8','O1','Oz','O2','F1','F2','Cz','AFz','POz'};

f=fields(variableStruct);

XLIMI = [-.100 .400];
yloglim = [10^-4  10^0];
    set(gcf,'Color','w')
    erppos = {[ 1 3 5],[2 4 6],[9 11 13],[10 12 14]} ;

    if strcmp(hipothesis,'faces') || strcmp(hipothesis,'NTF_NTO')
        chanplot = [15 9 10 12] %NTFvsNTO
    elseif strcmp(hipothesis,'VSNT_EX') || strcmp(hipothesis,'VSNTpre_VSNTpost') || strcmp(hipothesis,'VSNTpost_EX')
        chanplot = [13 14 10 12] % ASvsEX
    elseif strcmp(hipothesis,'target_AS')
        chanplot = [2 15 8 11];
    end
     %13 14 10 12  ]

for ifig = 1:4
    i = chanplot(ifig)
    er = subplot(8,2,erppos{ifig})
    % boxplot(h(i))
    ch = find(ismember({chanlocs.labels},chan{i}));
    
    for var = 1:length(variable2plot)
        struct(var).x = []
    end

%     x1 = [];
%     x2 = [];
    for iSuj = [1:16]
        %times   = variableStruct(iSuj).times;%(1:end-1);
        for var = 1:length(variable2plot)
            struct(var).x = [struct(var).x; variableStruct(iSuj).(variable2plot{var})(ch,:)];
        end
%         x1 = [x1; variableStruct(iSuj).(variable2plot{1})(ch,:)];
%         x2 = [x2; variableStruct(iSuj).(variable2plot{2})(ch,:)];
    end
    
    for var = 1:length(variable2plot)
        struct(var).X = squeeze(mean(struct(var).x,1));
        struct(var).eX = squeeze(std(struct(var).x,0,1))/sqrt(size(struct(var).x,1));
        struct(var).base = mean(struct(var).X(-.2<times & times<0));
        if baselineCorr==1
            struct(var).X = struct(var).X - struct(var).base
        end
        
    end   
    keyboard
%     X1 = squeeze(mean(x1,1));
%     X2 = squeeze(mean(x2,1));
%     eX1 = squeeze(std(x1,0,1))/sqrt(size(x1,1));
%     eX2 = squeeze(std(x2,0,1))/sqrt(size(x1,1));
 
%     base1 = mean(X1(-.2<times & times<0));
%     base2 = mean(X2(-.2<times & times<0));
%     if baselineCorr==1
%         X1    = X1-base1;
%         X2    = X2-base2;
%     end
    ylabel(strcat(chan{i},' [\muV]'),'FontSize',22) 

    box on                                    
    hold on
    %title(chan{i})
    
    
    for var = 1:length(variable2plot)
        struct(var).ylim = max([abs(max(abs(struct(var).X-struct(var).eX),[],'all')),abs(max(abs(struct(var).X+struct(var).eX),[],'all'))]);
    end
    %ylim = max([abs(max(abs(X1-eX1),[],'all')),abs(max(abs(X1+eX1),[],'all')),abs(max(abs(X2-eX2),[],'all')),abs(max(abs(X2+eX2),[],'all'))]);
    
    ylim_ = max([struct.ylim])
    YLIMI = [-ylim_-ylim_/2 ylim_];
    
    colors = {'k-' 'g-' 'r-' 'b-'};
    x = [times,     times(end:-1:1),            times(1)];
    for var = 1:length(variable2plot)
        y = [struct(var).X+struct(var).eX,    struct(var).X(end:-1:1)-struct(var).eX(end:-1:1), struct(var).X(1)+struct(var).eX(1)];
        patch(x,y,colors{var},'EdgeColor','none','FaceAlpha',0.5)
    end

%     % shade for se of cond1
%     x = [times,     times(end:-1:1),            times(1)];
%     y = [X1+eX1,    X1(end:-1:1)-eX1(end:-1:1), X1(1)+eX1(1)];
%     patch(x,y,'b','EdgeColor','none','FaceAlpha',0.5)
% 
%     % shade for se of cond2
%     x = [times,     times(end:-1:1),            times(1)];
%     y = [X2+eX2,    X2(end:-1:1)-eX2(end:-1:1), X2(1)+eX2(1)];
%     patch(x,y,'r','EdgeColor','none','FaceAlpha',0.5)

    % shade for time windows of interest
    x = [ventanas{1} ventanas{1}(end:-1:1) ventanas{1}(1) ];
    y = [8 8 -8 -8 8];
    %patch(x,y,'black','EdgeColor','none','FaceAlpha',0.4)

    x = [ventanas{2} ventanas{2}(end:-1:1) ventanas{2}(1) ];
    y = [8 8 -8 -8 8];
    %patch(x,y,'black','EdgeColor','none','FaceAlpha',0.4)

    x = [ventanas{3} ventanas{3}(end:-1:1) ventanas{3}(1) ];
    y = [8 8 -8 -8 8];
    %patch(x,y,'black','EdgeColor','none','FaceAlpha',0.4)

    x = [ventanas{4} ventanas{4}(end:-1:1) ventanas{4}(1) ];
    y = [8 8 -8 -8 8];
    %patch(x,y,'black','EdgeColor','none','FaceAlpha',0.4)
    
    h = [];
    hlegend = [];
    for var = 1:length(variable2plot)
        h(var) = plot(times,struct(var).X,colors{var},'LineWidth',3,'DisplayName',variable2plot{var});
        hlegend = [hlegend h(var)];
        if ~strcmp(variable2plot{var},'Intercept')
            if sum(strcmp(variable2plot,'Intercept'))>0
                p_block = find([tfce{var-1,1}.P_Values(ch,:)]<0.05);
            else
                p_block = find([tfce{var,1}.P_Values(ch,:)]<0.05);
            end
            difp = diff(p_block);
            indcut=1;
            for j = 1:length(difp)
                if difp(j) ~= 1
                    indcut = [indcut j];
                end
            end
            indcut = [indcut length(difp)];
            if sum(strcmp(variable2plot,'Intercept'))>0
                Y = tfce{var-1,1}.P_Values(ch,:); Y(Y>0.05)=1;
            else
                Y = tfce{var,1}.P_Values(ch,:); Y(Y>0.05)=1;
            end
            % imagesc(times(-99<times & times<1000*XLIMI(2)-1),[-0.97*var*ylim -0.9*var*ylim],Y(find(-99<times & times<1000*XLIMI(2)-1)),[0 0.05]); colormap hot;
            imagesc(times(-99<times & times<1000*XLIMI(2)-1),[-5 -4.5],Y(find(-99<times & times<1000*XLIMI(2)-1)),[0 0.05]); colormap hot;
        end
    end
%     h(1) = plot(times,X1,'b-','LineWidth',3);
%     h(2) = plot(times,X2,'r-','LineWidth',3,'DisplayName',variable2plot{2});
    
     if ifig == 1
        %legend(hlegend);
         legend([h(1) h(2)],variable2plot{1},variable2plot{2},'','')
     end

    plot([0 0],[-6,6],'k-')
    plot(XLIMI,[0 0],'k-')
    ax=gca;
    ax.FontSize = 18; 

    hold off
    
    % set(gca,'YLim',1.5*YLIMI,'XLim',XLIMI,'fontsize',15)
    currentPos = get(er,'Position');
    set(er,'Position',currentPos - [0,0,0,0.055])

%     str = sprintf('%s with %d  fixations and %s with %d fixations, 16 subjects', [variable2plot{1}(1:3)]...
%     ,sum([variableStruct.(['n_' cond{1}])]), [variable2plot{2}(1:3)],sum([variableStruct.(['n_' cond{2}])]));
%     str3=sprintf([model2analyse ' and fix >0.1'])

    a = get(gca,'XTickLabel');
    if ifig > 2
        xlabel('time [s]','FontSize',18);
    end
    
    xlim([-0.1,0.4])
    ylim([-6,6])
    set(gca, 'PlotBoxAspectRatio', [1,1,1])
end     

%% plot topoplots
variable2plot =  {'NCRank' 'NCRank_VSNT_EX'} %{'faces' hipothesis}% {'faces' hipothesis}%{'Intercept' hipothesis}% {hipothesis 'NCRank'} %{'Intercept' 'faces'}%    {'StandFixRank_VSNT_EX' 'StandFixRank_isEX'} % %{'Intercept' 'faces'}
ventanas={[0.05 0.15],[0.15 0.25],[0.25 0.35]}; %   %   ventanas={[0.05 0.1],[0.1 0.15],[0.14 0.25],[0.25 0.35]};  
    
% mdf 2/3/21 para señalar el índice de filas del gráfico de TOPOs sobre
% condiciones
jotas = [0 length(ventanas)] ;
% mdf 2/3/21 for loop sobre las ventanas de tiempos
[datavector1,datavector2] = deal([]);

if baselineCorr==1
    ventana_baseline = find(times>-0.2 & times<0);
    base1 = [];
    base2 = [];
    for iSuj = 1:Nsuj
        base1 = [base1, mean(variableStruct(iSuj).(variable2plot{1})(:,ventana_baseline),2)];
        base2 = [base2, mean(variableStruct(iSuj).(variable2plot{2})(:,ventana_baseline),2)];
    end
end

for i=1:length(ventanas)
    x1 = [];
    x2 = [];
    ventana = find(times>ventanas{i}(1) & times<ventanas{i}(2));

    min_max = min(length(variableStruct(1).(variable2plot{1})),length(variableStruct(1).(variable2plot{2})));
    ventana = ventana(ventana<=min_max);

    for iSuj = 1:Nsuj
        x1 = [x1, mean(variableStruct(iSuj).(variable2plot{1})(:,ventana),2)];
        x2 = [x2, mean(variableStruct(iSuj).(variable2plot{2})(:,ventana),2)];

    end
    
    if baselineCorr==1
        x1 = x1 - base1;
        x2 = x2 - base2;
    end    
    

    for j=jotas  
        if j==0
            datavector1 = [datavector1 mean(x1,2)]; % intercept 
        elseif j==length(ventanas)      
            datavector2 = [datavector2 mean(x2,2)];  % hipothesis
        end
    end
end
        
% max_pot1 = max(datavector1,[],'all');
% min_pot1 = min(datavector1,[],'all');
% abs_max1 = max(max_pot1, abs(min_pot1));
% 
% max_pot2 = max(datavector2,[],'all');
% min_pot2 = min(datavector2,[],'all');
% abs_max2 = max(max_pot2, abs(min_pot2));
figure(101); clf
[ha, pos] = fp.tight_subplot(2,length(ventanas),[0.0001 .1],[.00001 .0001],[0.1 0.1]); %[ha, pos] = fp.tight_subplot(2,length(ventanas),[.01 .1],[.01 .1],[0.1 0.1]);
for i=1:length(ventanas)
    for j=jotas 
        hold on 
        
        axes(ha(j+i));
        
        if i==1
            set(gca,'Visible','off');
            po=text(.5,.5,variable2plot{find(jotas==j)});
            set(po,'Position',[-1.7,0,0],'FontSize',21);
        end              
        if freqFilteredEEG        
            % topoplot(datavector1(:,[i]),chanlocs,'maplimits',[min_pot1,max_pot1],'electrodes','off');  colormap('gray');
            if contains(variable2plot{find(jotas==j)},'ank')
                extremos = [-6,6];
            else
                extremos = [-3,3];
            end
        else
            % topoplot(datavector1(:,[i]),chanlocs,'maplimits',[-abs_max1,abs_max1],'electrodes','off'); colormap(redblue(255))
            if contains(variable2plot{find(jotas==j)},'ank')
                extremos = [-6,6];
            else
                extremos = [-3,3];
            end
        end 
        if j==0 
             str1 = sprintf('%d' ,1000*ventanas{i}(1));
             str2 = sprintf('%d ms' ,1000*ventanas{i}(2));
             title(ha(i),[str1 '-' str2],'FontSize',30);%[ventanas{i}(1) ' ms' ;ventanas{i}(2) ' ms'] )
             topoplot(datavector1(:,[i]),chanlocs,'maplimits',extremos,'electrodes','off'); colormap(redblue(255))
        else
            topoplot(datavector2(:,[i]),chanlocs,'maplimits',extremos,'electrodes','off'); colormap(redblue(255))
        end
        
           
        
%         cb= colorbar;
%         %cb.Position = cb.Position + [.1+i*0.02, 0.1+j*0.1, 0.001, .04];
%         ylabel(cb, 'uV','FontSize',14);
%         set(gcf,'Color','w');
        if i==length(ventanas)
            % mdf 2/3/21 barra de referencia a la temperatura
            cb= colorbar;
            cb.Position = cb.Position + [.1, -0.1, 0.015, .2];
            cb.FontSize = cb.FontSize*4;
            ylabel(cb, 'uV','FontSize',30)
        end

        hold off
    end
end


% set(gcf,'Color','w')
%set(gcf,'Color','w','Position', get(0, 'Screensize'))
%saveas(gcf, '/media/cbclab/MARIADAFON1T/Analysis_2020/plots/eps/topoplots_hard&easy.eps','eps2c')

%% Bar plot number of trials per participant
figure(103); clf

Ncond = {strcat('n_',cond{1}), strcat('n_',cond{2})};
Ncond1 = [variableStruct.(['n_' cond{1}])];
Ncond2 = [variableStruct.(['n_' cond{2}])];

Ntrials = [Ncond1(:),Ncond2(:)];


bar(Ntrials);
yline(20);
legend(cond{1},cond{2});
xlabel('subject number','FontSize',18);
ylabel('amount of events','FontSize',18);

%% PO8 and PO7 average window 150-190 ms

window_time = times(0.149<=times & times<=0.191);
window_time2 = times(-0.1992<=times & times<=0.3984);

PO7_ampl = [];
PO8_ampl = [];
for iSuj = [1:16]
    PO7_ampl = [PO7_ampl; variableStruct(iSuj).NTF_NTO(25,find(ismember(times,window_time)))];
    PO8_ampl = [PO8_ampl; variableStruct(iSuj).NTF_NTO(62,find(ismember(times,window_time)))];
end

PO7_mean = mean(PO7_ampl,2);
PO8_mean = mean(PO8_ampl,2);

mean_str.PO7 = PO7_mean;
mean_str.PO8 = PO8_mean;

%% frontal and occipital channels average window 250-350 ms

iniTime = 0.15;
finTime = 0.25;

window_time = times(iniTime<=times & times<=finTime);
window_time_baseline = times(-0.1992<=times & times<0);

Fz_ampl = [];
Oz_ampl = [];
Fz_base = [];
Oz_base = [];
for iSuj = [1:16]
    Fz_ampl = [Fz_ampl; variableStruct(iSuj).NCRank(38,find(ismember(times,window_time)))];
    Oz_ampl = [Oz_ampl; variableStruct(iSuj).NCRank(29,find(ismember(times,window_time)))];
    if baselineCorr==1
        Fz_base = [Fz_base; variableStruct(iSuj).NCRank(38,find(ismember(times,window_time_baseline)))];
        Oz_base = [Oz_base; variableStruct(iSuj).NCRank(29,find(ismember(times,window_time_baseline)))];
    end
end

Fz_mean = mean(Fz_ampl,2);
Oz_mean = mean(Oz_ampl,2);
if baselineCorr==1
    Fz_mean_base = mean(Fz_base,2);
    Oz_mean_base = mean(Oz_base,2);
    Fz_mean = Fz_mean - Fz_mean_base;
    Oz_mean = Oz_mean - Oz_mean_base;
end

mean_str.Fz = Fz_mean;
mean_str.Oz = Oz_mean;

IT = int2str(iniTime*1000);
FT = int2str(finTime*1000);
fileOut = [folderOut,  '/FzOz_rank', IT, '_', FT, '.mat'];
save(fileOut,'mean_str')