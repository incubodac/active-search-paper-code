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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Edit this lines to choose a model to fit It should be written in Wilkinson
%notation, categorical and variables model with splines should be indicated
%cat(categorical variable), spl(variable,n)

%jek 11/11/2022: The configuration of the models is done in a separate
%file, because there are several parameters and comments
kk=1;
% figure(10); clf; set(gcf,'Color','w')
for RR=5:15
    jk2022_config_models;

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
    
    variable2TFCE   = 'NCRank_VSNT_EX';
    TFCEfilename = ['tfceResults_' variable2TFCE '.mat' ];
    fileTFCE  = fullfile(folderOut, TFCEfilename);
    load(fileTFCE)
    
    YY(kk)=sum(sum(tfceResult.P_Values(:,times>0)<0.05));kk=kk+1;
    
%     Y = tfceResult.P_Values;
%     figure(10 + ii); 
%         imagesc(times,1:64,Y([1:32 64:-1:33],:),[0 0.05]); 
%         
%         % Colormap
%         colormap hot; 
%         c=colorbar; 
%         ylabel(c,'p value','fontsize',42); 
%         set(c,'fontsize',35)
%     
%         xlabel('Time [s]','fontsize',42)
%         ylabel('Channels','fontsize',42)
%         title(variable2TFCE,'fontsize',42)
%         yticks([7 20 27 33 40 55])
%         nuChans = {chanlocs([1:32 64:-1:33]).labels};max
%         yticklabels(nuChans([7 25 27 33 40 62]))
%         set(gcf,'Color','w')
%         set(gca,'fontsize',35)
%         xlim([-0.1,0.4])
%         xline(xmin_sig,'Color',[0 0 0.5]+0.005,'linewidth', 5)
%         xline(0,'-k','linewidth', 5)
%         xline(xmax_sig,'Color',[0 0 0.5]+0.005,'linewidth', 5)
%         ax = gca;
%         get(gca)
%         set(gca,'linewidth', 5)
end

%%
NN = numel(tfceResult.P_Values);
figure(1); clf
    set(gcf,'Color','w')
    subplot(2,1,1)
        plot(5:15, log10(YY+1) ,'k.-','linewidth',3)
        set(gca,'YTick',[log10(1),log10(10),log10(100),log10(1000)],'YTickLabel',[1,10,100,1000])
        set(gca,'FontSize',10)
        xlabel('max rank')
        ylabel('sum significant samples (log10)')

    subplot(2,1,2)
        plot(5:15, YY ,'k.-','linewidth',3)
        set(gca,'FontSize',10)
        xlabel('max rank')
        ylabel('sum significant samples')


