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
% init_unfold

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

chanlocs = load('locationFile.mat');
chanlocs = chanlocs.e_loc;

% limits of significant time-window defined by visual inspection
xmin_sig = 0.0;
xmax_sig = 0.0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Edit this lines to choose a model to fit It should be written in Wilkinson
%notation, categorical and variables model with splines should be indicated
%cat(categorical variable), spl(variable,n)

%jek 11/11/2022: The configuration of the models is done in a separate
%file, because there are several parameters and comments
kk=1;
figure(10); clf; set(gcf,'Color','w','Position',[73 1 1848 961]); nn=0;
FS = 15; LW = 2;
me = 0.06; sp = 0.168; mi = 0.010;
% left bottom width heigth
pos = []; for jj = 5:-1:1; for ii = 1:5; pos = [pos; [me + (ii-1)*(sp+mi), me + (jj-1)*(sp+mi),sp,sp]]; end; end
pltchans = [7 20 27 33 40 55];
% variblesnames = {'(Intercept)','Task','Category','Rank','Task:Rank'};
variblesnames = {'Category','Task','Rank','Task:Rank'};


rr = 0;
for RR=5:15
    jk2022_config_models;

    if ismember(RR,[8 9 10 11 12])
        rr = rr+1;
        variables       = fieldnames(variableStruct);
        variables(cellfun(@(x) strcmp(x,'times'),variables)) = [];
        variables(cellfun(@(x) strcmp(x,'Suj'),variables)) = [];
%         variables = {'VSNT_EX','faces','NCRank','NCRank_VSNT_EX'};

        for ii = 1:length(variables)
            variable2TFCE   = variables{ii};
            [variableStruct, tfceResult] = load_data(who_is_running,session_path,datasetID, ...
                freqFilteredEEG,currentBand,modelName, ...
            variable2TFCE);        
            times   = variableStruct.times; 
        
            Y = tfceResult.P_Values;
%             nn=nn+1;
            nn = rr + 5*(ii-1);
            subplot(5,5,nn);
                imagesc(times,1:64,Y([1:32 64:-1:33],:),[0 0.05]); 
                
                % Colormap
                colormap hot; 
%                 if mod(nn,5)==0
                if nn==25
                    c=colorbar; 
                    ylabel(c,'p value','fontsize',FS); 
                    set(c,'fontsize',FS)
                end

%                 if mod(nn,5)==1
                if nn==21
                    ylabel('Channels','fontsize',FS)
                    nuChans = {chanlocs([1:32 64:-1:33]).labels};
                    yticks(pltchans)
                    yticklabels(nuChans(pltchans))
                else
                    yticks(pltchans)
                    yticklabels([])
                end
                if ismember(nn,1:5)
                    title(sprintf('max rank = %d',RR),'fontsize',FS)
%                     title(variblesnames{ii},'fontsize',FS)
                end
                if ismember(nn,21:25)
                    xlabel('Time [s]','fontsize',FS)
                    xticks(0:.100:.300)
                else
                    xticks(0:.100:.300)
                    xticks([])                    
                end
                set(gcf,'Color','w')
                set(gca,'fontsize',FS)
                xlim([-0.1,0.4])
%                 xline(xmin_sig,'Color',[0 0 0.5]+0.005,'linewidth', LW)
%                 xline(xmax_sig,'Color',[0 0 0.5]+0.005,'linewidth', LW)
                xline(0,'-k','linewidth', LW)
                ax = gca;
%                 get(gca)
                set(gca,'Position', pos(nn,:))
                set(gca,'linewidth', LW)

            if (ii==5 && RR==8);    Y8int   = sum(Y(:,times>.1)<0.05,2); Y8int = Y8int/max(Y8int); end
            if (ii==5 && RR==10);   Y10int  = sum(Y(:,times>.1)<0.05,2); Y10int = Y10int/max(Y10int); end
            if (ii==5 && RR==12);   Y12int  = sum(Y(:,times>.1)<0.05,2); Y12int = Y12int/max(Y12int); end
        end
    end

    variable2TFCE   = 'NCRank_VSNT_EX';
    [variableStruct, tfceResult] = load_data(who_is_running,session_path,datasetID, ...
                    freqFilteredEEG,currentBand,modelName, ...
                    variable2TFCE);        
    times   = variableStruct.times; 
    YY(kk)=sum(sum(tfceResult.P_Values(:,times>0)<0.05));
    kk=kk+1;
end

%%


path_data = '/home/juank/repos/corregistro/Analysis_dac/data/input/';
path2save = '/home/juank/Desktop/EJN2022/data/output/';

df = fun_importfile_EJ2022_cvs([path2save,'EJN2022_16subj_GeneralFilters.csv']);
df_absent = fun_importfile_EJ2022_cvs([path2save,'EJN2022_16subj_VSNTabsent_EX_manyNranks.csv']); 
df_pretar = fun_importfile_EJ2022_cvs([path2save,'EJN2022_16subj_VSNTpre_EX_manyNranks.csv']);

Xlabels = {'EX','VSabsent','VSpretarget'};

% gcapos = [.1 .2 .525 .6; .65 .525 .25 .275; .65 .2 .25 .275]
gcapos = [.1 .525 .25 .275; .1 .2 .25 .275; .375 .2 .525 .6];
NN = numel(tfceResult.P_Values);
figure(20); clf
    set(gcf,'Color','w','Position',[73         625        1847         337]); 
    subplot(2,4,5)
        hx = 1.5:20.5;
        cols = [.75,.1,.1;.1,.75,.1;.1,.1,.75];
        hold on
            X = {rank_trials(df(df.trial_type=='VS',:)), ...
                    rank_trials(df_absent(df_absent.trial_type=='VS',:)),...
                    rank_trials(df_pretar(df_pretar.trial_type=='VS',:))};
            for ii = 1:3
                fprintf('%s: len = %d, median = %d, max = %d\n',Xlabels{ii},length(X{ii}),median(X{ii}),max(X{ii}));
                hy = hist(X{ii}, hx);%, histtype='step', stacked=True, fill=False, label=Xlabels[i], density=True)
                plot(hx,hy/sum(hy),'-',color=cols(ii,:),linewidth=LW);
            end
        hold off
        box on
        set(gca,'XAxisLocation','bottom')
    %     set(gca,'YAxisLocation','right')
        set(gca,'xlim',[0 20])
        set(gca,'FontSize',FS)
        xlabel('trial length')
        ylabel('Prob.')
    %         legend({'EX','VS-Absent','VS-Pretarget'}); legend boxoff
        set(gca,'linewidth', LW)
        set(gca,'Position',gcapos(2,:))

    subplot(2,4,1)
        hx = 1.5:20.5;
        cols = [.75,.1,.1;.1,.75,.1;.1,.1,.75];
        hold on
            X = {   df(df.trial_type=='VS',:).rank, ...
                    df_absent(df_absent.trial_type=='VS',:).rank,...
                    df_pretar(df_pretar.trial_type=='VS',:).rank};
            for ii = 1:3
                fprintf('%s: len = %d, median = %d, max = %d\n',Xlabels{ii},length(X{ii}),median(X{ii}),max(X{ii}));
                hy = hist(X{ii}, hx);%, histtype='step', stacked=True, fill=False, label=Xlabels[i], density=True)
                plot(hx,hy/sum(hy),'-',color=cols(ii,:),linewidth=LW);
            end
        hold off
        box on
        set(gca,'XAxisLocation','top')
    %     set(gca,'YAxisLocation','right')
        set(gca,'xlim',[0 20])
        set(gca,'FontSize',FS)
        xlabel('fixation rank')
        ylabel('Prob.')
        legend({'EX','VS-Absent','VS-Pretarget'}); legend boxoff
        set(gca,'linewidth', LW)
        set(gca,'Position',gcapos(1,:))

    subplot(2,4,8)
        hold on
            xline(8,'--k')%,'linewidth', LW)
            xline(9,'--k')%,'linewidth', LW)
            xline(10,'--k')%,'linewidth', LW)
            xline(11,'--k')%,'linewidth', LW)
            xline(12,'--k')%,'linewidth', LW)
            plot(5:15, YY ,'k.-','linewidth', 5)
        hold off
        box on
        set(gca,'XAxisLocation','top')
        set(gca,'YAxisLocation','right')
        set(gca,'xlim',[4 16])
        set(gca,'FontSize',FS)
        xlabel('max. fixation rank')
        ylabel('TRF signif samples (10^3)')
        set(gca,'ytick',[1000 3000],'yticklabel',{'1.0' '3.0'})
        set(gca,'linewidth', LW)
        set(gca,'Position',gcapos(3,:))
        
%%
% path_data = '/home/juank/repos/corregistro/Analysis_dac/data/input/';
% path2save = '/home/juank/Desktop/EJN2022/data/output/';
% 
% df = fun_importfile_EJ2022_cvs([path2save,'EJN2022_16subj_GeneralFilters.csv']);
% df_absent = fun_importfile_EJ2022_cvs([path2save,'EJN2022_16subj_VSNTabsent_EX_manyNranks.csv']); 
% df_pretar = fun_importfile_EJ2022_cvs([path2save,'EJN2022_16subj_VSNTpre_EX_manyNranks.csv']);
% 
% Xlabels = {'EX','VSabsent','VSpretarget'};
% 
% NN = numel(tfceResult.P_Values);
% figure(20); clf
%     set(gcf,'Color','w','Position',[675 9 430 953])
%     subplot(8,3,4); topoplot(Y8int,chanlocs,'maplimits',[0 1]); colormap('gray'); %colorbar;
%     subplot(8,3,5); topoplot(Y10int,chanlocs,'maplimits',[0 1]); colormap('gray'); %colorbar;
%     subplot(8,3,6); topoplot(Y12int,chanlocs,'maplimits',[0 1]); colormap('gray'); %colorbar;
% 
%     subplot(4,1,2)
%         hold on
%             xline(8,'--k','linewidth', LW)
%             xline(10,'--k','linewidth', LW)
%             xline(12,'--k','linewidth', LW)
%             plot(5:15, YY ,'k.-','linewidth', LW)
%         hold off
%         box on
%         set(gca,'xlim',[0 20])
%         set(gca,'FontSize',FS)
%         xlabel('max rank')
%         ylabel('# signif samples (10^3)')
%         set(gca,'ytick',[1000 3000],'yticklabel',{'1.0' '3.0'})
% 
%     subplot(4,1,3)
%         hx = 1.5:20.5;
%         cols = [.75,.1,.1;.1,.75,.1;.1,.1,.75];
%         hold on
%             X = {rank_trials(df(df.trial_type=='VS',:)), ...
%                     rank_trials(df_absent(df_absent.trial_type=='VS',:)),...
%                     rank_trials(df_pretar(df_pretar.trial_type=='VS',:))};
%             for ii = 1:3
%                 fprintf('%s: len = %d, median = %d, max = %d\n',Xlabels{ii},length(X{ii}),median(X{ii}),max(X{ii}));
%                 hy = hist(X{ii}, hx);%, histtype='step', stacked=True, fill=False, label=Xlabels[i], density=True)
%                 plot(hx,hy/sum(hy),'-',color=cols(ii,:),linewidth=LW);
%             end
%         hold off
%         box on
%         set(gca,'xlim',[0 20])
%         set(gca,'FontSize',FS)
%         xlabel('trial length')
%         ylabel('P(length)')
% %         legend({'EX','VS-Absent','VS-Pretarget'}); legend boxoff
% 
%     subplot(4,1,4)
%         hx = 1.5:20.5;
%         cols = [.75,.1,.1;.1,.75,.1;.1,.1,.75];
%         hold on
%             X = {   df(df.trial_type=='VS',:).rank, ...
%                     df_absent(df_absent.trial_type=='VS',:).rank,...
%                     df_pretar(df_pretar.trial_type=='VS',:).rank};
%             for ii = 1:3
%                 fprintf('%s: len = %d, median = %d, max = %d\n',Xlabels{ii},length(X{ii}),median(X{ii}),max(X{ii}));
%                 hy = hist(X{ii}, hx);%, histtype='step', stacked=True, fill=False, label=Xlabels[i], density=True)
%                 plot(hx,hy/sum(hy),'-',color=cols(ii,:),linewidth=LW);
%             end
%         hold off
%         box on
%         set(gca,'xlim',[0 20])
%         set(gca,'FontSize',FS)
%         xlabel('trial length')
%         ylabel('P(length)')
%         legend({'EX','VS-Absent','VS-Pretarget'}); legend boxoff
        
%%
% NN = numel(tfceResult.P_Values);
% figure(21); clf
%     set(gcf,'Color','w')
%     subplot(2,1,1)
%         plot(5:15, log10(YY+1) ,'k.-','linewidth',LW)
%         set(gca,'YTick',[log10(1),log10(10),log10(100),log10(1000)],'YTickLabel',[1,10,100,1000])
%         set(gca,'FontSize',FS)
%         xlabel('max rank')
%         ylabel('sum significant samples (log10)')
% 
%     subplot(2,1,2)
%         plot(5:15, YY ,'k.-','linewidth',LW)
%         set(gca,'FontSize',FS)
%         xlabel('max rank')
%         ylabel('sum significant samples')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function max_rank = rank_trials(df)
    trial_subject = df.trial_subject;
    trial_number = df.trial_number;
    fixrank = df.rank;
    unique_sub = unique(trial_subject);
    max_rank = [];
    for indsub = 1:length(unique_sub)
        sfixrank        = fixrank(trial_subject==unique_sub(indsub));
        strial_number   = trial_number(trial_subject==unique_sub(indsub));
        unique_trial = unique(strial_number);
        for indtrial = 1:length(unique_trial)
            tsfixrank = sfixrank(strial_number==unique_trial(indtrial));
            max_rank = [max_rank, max(tsfixrank)];
        end
    end
end

function [variableStruct, tfceResult] = load_data(who_is_running,session_path,datasetID, ...
                    freqFilteredEEG,currentBand,modelName, ...
                    variable2TFCE)
    if strcmp(who_is_running,'M');      auxdir      =  session_path.fixs;           % fixs
    elseif strcmp(who_is_running,'D');  auxdir      =  session_path.aux_data;       % Analysis_aux_data
    elseif strcmp(who_is_running,'J');  auxdir      =  session_path.fixs;           % python: the csv with the output of select_fixs_unfold-EJN2022.ipynb
    end
    
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
    
    variables       = fieldnames(variableStruct);
    variables(cellfun(@(x) strcmp(x,'times'),variables)) = [];
    variables(cellfun(@(x) strcmp(x,'Suj'),variables)) = [];
    
    TFCEfilename = ['tfceResults_' variable2TFCE '.mat' ];
    fileTFCE  = fullfile(folderOut, TFCEfilename);
    load(fileTFCE)
    
end


