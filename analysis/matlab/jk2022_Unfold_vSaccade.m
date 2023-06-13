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
who_is_running = 'D'; % 'M': Maria or Damian or Juan % jek 10/11/22
sanity_checks = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(who_is_running,'M')
    cd  /home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/functions
elseif strcmp(who_is_running,'D')
    cd ~/Documents/repos/corregistro/EJN2022/matlab/    %   dac 1/5/23 directory
elseif strcmp(who_is_running,'J')
    cd /home/juank/Desktop/EJN2022/matlab/                              % jek 10/11/22
end

if strcmp(who_is_running,'M');      addpath(genpath('/home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/UNFOLD'))
elseif strcmp(who_is_running,'D');  addpath(genpath('/Users/dac/Documents/GitLab/corregistro/Analysis_dac/UNFOLD')) % mdf 8/6/21 DAC directory
elseif strcmp(who_is_running,'J');  addpath(genpath('/home/juank/repos/corregistro/Analysis_dac/UNFOLD'))           % jek 10/11/22
end

[runpath,code_path,session_path,whoisrunning]=add_paths_EEG_ET_beta(0,'D');% edited by jek 10/11/22
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
RR=5; jk2022_config_models;


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
    fprintf('Temporal EEG files: %s\n',folderPrepared)

% Models output folder for UNFOLD (remaining STEPs)        
models = fullfile(folderUnf, [datasetID, 'models/'] );
    if ~exist(models, 'dir'); mkdir(models); end

if freqFilteredEEG; folderModelName     = [currentBand.name '/' modelName '/'];
else;               folderModelName     = [modelName '/'];
end

folderOut               = fullfile(models,folderModelName);
    if ~exist(folderOut, 'dir'); mkdir(folderOut); end
    fprintf('Output model files: %s\n',folderOut)

fprintf('****************************************************************\n');
%% 1st STEP Create EEG with feature information embeded into EEGLAB struct
% This only needs to be run once (for each CSV)
if 1
    %%%%%%%%%%%%%%%%
    epochedByFixs = 0;
    useCsv        = 1; % jek 11/11/2022: Never used
    %%%%%%%%%%%%%%%%
    
    filename = [auxdir fixs_csv_file]; % csv file with fixations selected
        if exist(filename, 'file')
          % File exists.  Read it into a table.
          Tall = readtable(filename);
        else
          % File does not exist.  Warn the user.
          errorMessage = sprintf('Error: file not found:\n\n%s', filename);
          uiwait(errordlg(errorMessage));
          return;
        end
    
    %mdf 22/6/21 convert this column because matlab readtable interpret columns
    %with nan values as string
    %Tall.(hipothesis) = arrayfun(@(x) str2double(x), Tall.(hipothesis));
    %suj = fp.cfg.subjname ;
    
    for iSuj = subj_index    
        if freqFilteredEEG; sujName = names{iSuj}(1:3);
        else;               sujName = names{iSuj}(7:9);
        end
        
        rows = strcmp(Tall.trial_subject , sujName);    % mdf 8/6/21 select the rows of the current subject
        T = Tall(rows,:);                               % table of fixations for one subject  
        fp.cfg.inFile = [inFolder names{iSuj}];
        inFile = fp.cfg.inFile;
        
        % tload_start = tic;
        EEG = pop_loadset(inFile);
        % tload_end = toc(tload_start)
    
        if freqFilteredEEG && baselineCorr %%% power for EEG filtered by freq-band %%% mdf 8/4/22
            DATA = EEG.data;
            % power [of abs(Hilbert)]
            DATA = DATA.^2;
            times   = EEG.times;
            Ntrials = size(DATA,3);
            % trial by trial baseline correction
            base1 = mean(DATA(:,min(times)<times & times<0,:),2);
            base1 = squeeze(base1);
            invBase = (1./base1); 
              for i=1:64
                  for j=1:length(times)
                    aux = DATA(i,j,:);
                    lala = squeeze(invBase(i,:)).*transpose(squeeze(aux));
                    DATA(i,j,:)=reshape(lala, [1,1,Ntrials]);
                    DATA(i,j,:) = 10*log10(DATA(i,j,:));
                    EEG.data(i,j,:) = DATA(i,j,:);
                  end
              end
        end
        
        if epochedByFixs == 1
            EEG = eeg_epoch2continuous(EEG);
            EEG = pop_epoch( EEG, {  'fixation'  },[-0.2 0.8], 'newname', 'fixs and conditions epoch', 'epochinfo', 'yes');
            EEG = pop_rmbase( EEG, [-199.2188 0]); 
        end
        
        % tcont_start = tic;
        EEG = eeg_epoch2continuous(EEG);                            % mdf 8/6/21 concatenate trials
        % tcont_end = toc(tcont_start)
        
        % tanexo_start = tic;
        EEG = fp.renameEventsUnfoldSaccade(EEG, T ,hipothesis, 0.08);  % analysis made so far was at .28 last before argument istarget == 0 for distractors only
        % tanexo_end = toc(tanexo_start)
        
        % tcheck_start = tic;
        nbchan = EEG.nbchan;
        if  nbchan > 64
           EEG = pop_select(EEG,'nochannel',64+1:nbchan);
           EEG = eeg_checkset(EEG);                                 % mdf 8/6/21 check that after changes the number of channels didn't change
           fprintf('Unwanted channels erased\n');
        end
        % tcheck_end = toc(tcheck_start)
        
        outFile = [sujName  '_prepared_E_H_unfold'] ;               % mdf 9/6/21 overwrite this files
        fprintf('saving %s EEG in %s folder \n', outFile, folderPrepared);
    
        pop_saveset(EEG,[folderPrepared  outFile]);
    end
end

%% 2nd STEP run model and create ERP matrix to analyse
names = fp.loadSubjects(folderPrepared, 1, 1, 'set');
unfoldResults =  [];

suj = fp.cfg.subjname;

for iSuj = 1:length(names)
    tstart = tic;
    fprintf('Runing model for subject %d: %.3f\n',iSuj, tstart)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Load EEG data: subject %d: %.3f\n',iSuj, toc(tstart))    
    fp.cfg.inFile = [folderPrepared names{iSuj}];
    inFile              = fp.cfg.inFile;
    EEG = pop_loadset(inFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Select fixations based on conditions: subject %d: %.3f\n',iSuj, toc(tstart))
    if filter_faces % select only fixations to faces
        fprintf('\tFilter faces\n')
        event_faces = [];
        for ev = 1:numel([EEG.event]) 
            %keyboard
            if EEG.event(ev).faces == 1
                event_faces = [event_faces ev];
            end
        end
        EEG.event(event_faces) = [];
    end
    if filter_EX    % select only fixations in the exploration trials
        fprintf('\tFilter EX\n')
        event_EX = [];
        for ev = 1:numel([EEG.event]) 
            %keyboard
            if EEG.event(ev).VSNT_EX == 0
                event_EX = [event_EX ev];
            end
        end
        EEG.event(event_EX) = [];
    end
    if filter_AS    % select only fixations in the active search trials
        fprintf('\tFilter AS\n')
        event_AS = [];
        for ev = 1:numel([EEG.event]) 
            %keyboard
            if EEG.event(ev).VSNT_EX == 1
                event_AS = [event_AS ev];
            end
        end
        EEG.event(event_AS) = [];
    end
    if ~filter_faces && ~filter_EX && ~filter_AS;
        fprintf('\tNo filter is applied\n')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Defining the design matrix: subject %d: %.3f\n',iSuj, toc(tstart))
%     tend_load = toc(tstart_load)    
%     latency = [EEG.event.latency];
%     perm = randperm(length(latency));
%     PermLatency = latency(perm);
%     for il = 1:length(EEG.event)
%         EEG.event(il).OrigLatency = latency(il)
%         EEG.event(il).PermLatency = PermLatency(il)
%     end
%     tstart_unfold = tic;
    
    cfgDesign = [];
%     cfgDesign.eventtypes = {'fixation'}; % we model the fixation onset
%     cfgDesign.formula = model ;
    cfgDesign.eventtypes = {{'fixation'},{'saccade'}}; % we model the fixation onset
    cfgDesign.formula = {model, 'y~1'} ;
    
    %One needs to be careful to not overfit the data, regularization or cross validation can help here.
    EEG = uf_designmat(EEG,cfgDesign);
    
    % DELETE NANS mdf 29/6/21 
    cfgNueva = [];
    cfgNueva.method = 'drop';
    EEG = uf_imputeMissing(EEG,cfgNueva);
    
    % DETECT NOISY DATA
    winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',250); % before350posterversion % mdf 8/6/21 to detect sporius channels with amplitude over...
    winrej = bad_ET_ArtifactDetect(EEG,winrej); %DAC 27 sept 2021 detect bad_ET marks
    
    %TIME EXPAND
    cfgTimeexpand = [];
    cfgTimeexpand.timelimits = [-.2,.4];
    EEG = uf_timeexpandDesignmat(EEG,cfgTimeexpand);
    
    %REJECT NOISY DATA
    EEG = uf_continuousArtifactExclude(EEG,struct('winrej',winrej)); 
    % we find data using a simple threshold tool. More complex algorithms or manual cleaning is recommended!

    % for il = 1:length(EEG.event); EEG.event(il).latency = EEG.event(il).PermLatency; end   
    EEG.event(ismember({EEG.event.type},{'bad_ET','boundary'})) = []; %DAC mod 30 sept 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Fitting the model: subject %d: %.3f\n',iSuj, toc(tstart))
    EEG= uf_glmfit(EEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Saving results of the model: subject %d: %.3f\n',iSuj, toc(tstart))    
    suje = names{iSuj}; sujName = suje(1:3);
    
    n = struct(condi{1},sum([EEG.event.(hipothesis)]==1),condi{2},sum([EEG.event.(hipothesis)]==0));
    EEG.unfold.N = struct(condi{1}, n.(condi{1}), condi{2},n.(condi{2})); % unreliable
    
    unfoldResults.(sujName) = EEG.unfold;
    
    if sanity_checks
        % folder_sanity = [folderOut 'afterBadETafterBoundary_onlyFix/'];
        if strcmp(hipothesis,'NTF_NTO')
            ch1 = 'PO7'
        elseif strcmp(hipothesis,'VSNT_EX')
            ch1 = 'F1'
        end
        fp.plotUnfoldErpimageCheck(EEG,ch1,hipothesis,'fixation',10)
        % saveas(gcf,'Barchart.png')
        if strcmp(hipothesis,'NTF_NTO')
            ch2 = 'PO8'
        elseif strcmp(hipothesis,'VSNT_EX')
            ch2 = 'F2'
        end
        fp.plotUnfoldErpimageCheck(EEG,ch2,hipothesis,'fixation',10)       
        fp.plotUnfoldErpimageCheck(EEG,'O1',hipothesis,'fixation',10)     
        fp.plotUnfoldErpimageCheck(EEG,'O2',hipothesis,'fixation',10)     
    end

    fileOut = [folderOut  'Results.mat'];
    save(fileOut,'unfoldResults')
    
end

%% 3rd STEP Generate Data matrix for TFCE and plots for variables of interest
if exist(folderOut, 'dir')
    load(fullfile(folderOut, 'Results.mat'));
    variableStruct =[];
    subjects   = fieldnames(unfoldResults);
    for iSuj  = 1:length(fieldnames(unfoldResults))
        suj   = subjects{iSuj}; 
        result = unfoldResults.(suj);
        variableStruct  = fun_makeStructToERPunfold(variableStruct,variablesOI,iSuj,suj,result);
    end
    fileOut   = fullfile(folderOut,'variableStruct.mat');
    save(fileOut,'variableStruct')
else
    fprintf("the model doesn't  exist, check model name or fit it \n")
end

%% 4th STEP run TFCE
fileStruct      = fullfile(folderOut,'variableStruct.mat');
load(fileStruct)
variables       = fieldnames(variableStruct);
variables(cellfun(@(x) strcmp(x,'times'),variables)) = [];
variables(cellfun(@(x) strcmp(x,'Suj'),variables)) = [];

% If you want to run a single variable
% variables       = {'VSNT_EX_faces'};

for ii = 1:length(variables)
    variable2TFCE   = variables{ii};
    if exist(folderOut, 'dir')
        load(fileStruct);
        times   = variableStruct.times; 
    %     chanlocs = load(fullfile(session_path.aux_data,'locationFile.mat'));
        chanlocs = load('locationFile.mat'); % jek 10/11/2022
        chanlocs = chanlocs.e_loc;
        if baselineCorr == 1 && ~freqFilteredEEG
            variableStruct  = fun_applyBaselineCorrection(variableStruct,variable2TFCE,freqFilteredEEG);
            fileOut         = fullfile(folderOut,['beta_', variable2TFCE,'_BC.mat']);
    
            save(fileOut,'variableStruct');    
        end
        tfceResult  = fun_runTFCE(fp,folderOut,variableStruct, variable2TFCE,chanlocs);
    else
        fprintf("the model doesn't  exist, check model name and variableStruct \n")
    end
end

%% Functions 

function ERPstruct  = fun_makeStructToERPunfold(ERPstruct,variablesOI,iSuj,suj,result)
    
    ERPstruct(iSuj).times    = result.times;
    ERPstruct(iSuj).Suj      = suj;
%     if sum(contains(variablesOI, hipothesis))>0
%         condi = split(hipothesis,'_');
%         ERPstruct(iSuj).(['n_' condi{1}]) = result.N.(condi{1});
%         ERPstruct(iSuj).(['n_' condi{2}]) = result.N.(condi{2});
%     end

    for var = 1:length(variablesOI)
        covariate             = find(ismember(result.variablenames,variablesOI{var}));
        if strcmp(variablesOI{var},'(Intercept)')
            ERPstruct(iSuj).Intercept = result.beta_dc(:,:,covariate);    
        else
            if contains(variablesOI{var},':')
                lala = replace(variablesOI{var},':','_');
                ERPstruct(iSuj).(lala) = result.beta_dc(:,:,covariate);
            else
                ERPstruct(iSuj).(variablesOI{var}) = result.beta_dc(:,:,covariate);      
            end
        end  
    end
end
        
function tfceResult  = fun_runTFCE(fp,folder,variableStruct,variable,chanlocs)
    TFCEvars = fieldnames(variableStruct);
    data1 = [];
    
    if strcmp(variable,'(Intercept)')
        variable = 'Intercept';
    end
    var = find(ismember(TFCEvars,variable));
    for iSuj = 1:length(variableStruct)
       data = variableStruct(iSuj).(TFCEvars{var});
       data1= cat(3,data1,data);
    end
    data1 = permute(data1,[3 1 2]);
    Summary = data1;

    % save([folder '/data1'],'Summary')
    data2 = zeros(size(data1));

    tfceResult = fp.ept_TFCE(data1,data2,chanlocs,'nPerm',2000, 'type','d', 'rSample',512);
    %tfceResult = fp.ept_TFCE()
    
    save([folder '/tfceResults_' variable],  'tfceResult')
end

function variableStruct = fun_applyBaselineCorrection(variableStruct,variable,freqFilteredEEG)
    
    TFCEvars = fieldnames(variableStruct);    
    if strcmp(variable,'(Intercept)')
        variable = 'Intercept';
    end
    var = find(ismember(TFCEvars,variable));                
    for iSuj = 1:numel(variableStruct)
          times   = variableStruct(iSuj).times;

          x1 = variableStruct(iSuj).(TFCEvars{var});

          base1 = mean(x1(:,min(times)<times & times<0),2);
          
          x1    = x1-base1;      
          
          variableStruct(iSuj).(TFCEvars{var})=  x1;           
    end
end
