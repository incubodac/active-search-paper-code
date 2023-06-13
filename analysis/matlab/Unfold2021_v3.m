%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/6/21 DAC                 %   
% 9/6/21 modified by MDF     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
who_is_running = 'M'; % 'M': Maria or Damian or Juan % jek 10/11/22
sanity_checks = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(who_is_running,'M')
    cd  /home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/functions
elseif strcmp(who_is_running,'D')
    cd /Users/dac/Documents/GitLab/corregistro/Analysis_dac/functions   % mdf 8/6/21 DAC directory
elseif strcmp(who_is_running,'J')
    cd /home/juank/Desktop/EJN2022/matlab/                              % jek 10/11/22
end

[runpath,code_path,session_path,whoisrunning]=add_paths_EEG_ET_beta(1,2);% edited by jek 10/11/22

eeglabpath = code_path.eeglabpath


addpath(eeglabpath); 

ftpath      =        code_path.toolbox ;
%addpath('/media/dac/data/repos/visualsearch_eegem/codes')
addpath(ftpath); 
addpath(genpath(code_path.corregistro_fn))
addpath(code_path.my_functions)
addpath(code_path.psyco )
addpath(runpath)
addpath(code_path.functions)

clear fp; fp = FP_epochAnalysis();

eeglab
matdir              = session_path.matfiles;
datapath            = session_path.data_analysis; 
% auxdir              =  session_path.aux_data; % mdf 8/6/21 DAC directory

cd(datapath)

if strcmp(who_is_running,'M')
    addpath(genpath('/home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/UNFOLD'))
elseif strcmp(who_is_running,'D')
    addpath(genpath('/Users/dac/Documents/GitLab/corregistro/Analysis_dac/UNFOLD')) % mdf 8/6/21 DAC directory
elseif strcmp(who_is_running,'J')
    addpath(genpath('/home/juank/repos/corregistro/Analysis_dac/UNFOLD'))           % jek 10/11/22
end

%
fp.cfg.subjname         =  session_path.subjname;
fp.cfg.sessionfilenames = session_path.sessionfilenames;
fp.cfg.matdir           = matdir;
folderUnf               = fullfile(session_path.out, 'Unfold2021/' ); 
folderPrepared          = fullfile(folderUnf, 'preparedEEG/' ); 

init_unfold

%Edit this lines to choose a model to fit It should be written in Wilkinson
%notation, categorical and variables model with splines should be indicated
%cat(categorical variable), spl(variable,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hipothesis = 'VSNT_EX';
modelName = 'Inter+task+faces+NCrank+interTaskNCRank';
model     = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank

%e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp' you can
variablesOI   = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'}; %'StandFixRank:isEX'};
%e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
cond = {'EX','AS'};%%%cond = split(hipothesis,'_');
fixs_csv_file = '16subj_VSNTpostvsEX_for_unfold_GlobalStandarisedCenteredRank.csv'%'16subj_hardVSeasy_for_unfold.csv';

baselineCorr    = 1;
% ourPermTest     = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose the freq-band
Theta = struct('limits',[4  8],'name','Theta');
Alpha = struct('limits',[8 13],'name','Alpha');
Beta = struct('limits',[15 25],'name','Beta');
LowBeta = struct('limits',[14 19],'name','LowBeta');
HighBeta = struct('limits',[20 30],'name','HighBeta');
LowGamma  = struct('limits',[30 40],'name','LowGamma');

currentBand = Theta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqFilteredEEG = 1; % 1 if we used the EEG filtered in freq
filter_faces = 0; % 1 if filter fix to faces O otherwise
filter_EX = 0;
filter_AS = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(who_is_running,'M')
    auxdir              =  session_path.fixs;
else
    auxdir              =  session_path.aux_data;
end

if freqFilteredEEG
    inFolder            =  [session_path.freqEEG currentBand.name '/'];
    names = fp.loadSubjects(inFolder, 1, 1, 'set');
    subj_index = [1:length(names)];
else   
    inFolder            =  session_path.data_analysis;
    names = fp.loadSubjects(inFolder, 1, 1, 'set');
    subj_index = [1:4 6:11 13:15 19 22:23];
end

Nsuj = length(subj_index);

%% 1st STEP Create EEG with feature information embeded into EEGLAB struct
%%%%%%%%%%%%%%%%
epochedByFixs = 0;
useCsv        = 1;
%%%%%%%%%%%%%%%%

if ~exist(folderUnf, 'dir')
    mkdir(folderUnf)    
end
if ~exist(folderPrepared, 'dir')
    mkdir(folderPrepared)    
end

if freqFilteredEEG
    folderOut                =  [folderPrepared currentBand.name '/'];
else 
    folderOut                = folderPrepared;
end

filename = [auxdir fixs_csv_file] % csv file with fixations selected

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
    
    if freqFilteredEEG
        sujName = names{iSuj}(1:3);
    else
        sujName = names{iSuj}(7:9);
    end
    
    rows = strcmp(Tall.trial_subject , sujName); % mdf 8/6/21 select the rows of the current subject

    T = Tall(rows,:); %table of fixations for one subject  
 
    fp.cfg.inFile = [inFolder names{iSuj}];

    inFile              = fp.cfg.inFile;
    
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
    % keyboard
    if epochedByFixs == 1
        EEG = eeg_epoch2continuous(EEG)
        EEG = pop_epoch( EEG, {  'fixation'  },[-0.2 0.8], 'newname', 'fixs and conditions epoch', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-199.2188 0]); 
        %keyboard
    end
    
    % tcont_start = tic;
    EEG = eeg_epoch2continuous(EEG) % mdf 8/6/21 concatenate trials
    % tcont_end = toc(tcont_start)
    
    % tanexo_start = tic;
    EEG = fp.renameEventsUnfold2021(EEG, T ,hipothesis, 0.08);%analysis made so far was at .28 last before argument istarget == 0 for distractors only
    % tanexo_end = toc(tanexo_start)
    
    %keyboard
    % tcheck_start = tic;
    nbchan = EEG.nbchan;
    if  nbchan > 64
       EEG = pop_select(EEG,'nochannel',64+1:nbchan);
       EEG = eeg_checkset(EEG); % mdf 8/6/21 check that after changes the number of channels didn't change
       fprintf('Unwanted channels erased\n')
    end
    % tcheck_end = toc(tcheck_start)
    
    outFile = [sujName  '_prepared_E_H_unfold'] ; % mdf 9/6/21 overwrite this files
    fprintf('saving %s EEG in %s folder \n', outFile, folderOut)

    pop_saveset(EEG,[folderOut  outFile ])
end

%% 2nd STEP run model and create ERP matrix to analyse

inFolder                =  fullfile(session_path.out, 'Unfold2021/' );
fp.cfg.subjname         =  session_path.subjname;
fp.cfg.sessionfilenames = session_path.sessionfilenames;
fp.cfg.matdir           = matdir;

if freqFilteredEEG
    folderPrepared          = fullfile(inFolder, 'preparedEEG/', currentBand.name, '/');
else   
    folderPrepared          = fullfile(inFolder, 'preparedEEG/' );
end

models = fullfile(inFolder, 'models/' );
if ~exist(models, 'dir')
    mkdir(models)
end

if freqFilteredEEG
    folderModelName         = [currentBand.name '/' modelName '/'];
else   
    folderModelName         = [ modelName '/'];
end

folderOut               = fullfile(models,folderModelName);

if ~exist(folderOut, 'dir')
    mkdir(folderOut)
end

names = fp.loadSubjects(folderPrepared, 1, 1, 'set');
unfoldResults =  [];

suj = fp.cfg.subjname ;

for iSuj = 1:length(names)
    iSuj
    % tstart_load = tic;
    
    fp.cfg.inFile = [folderPrepared names{iSuj}];
    inFile              = fp.cfg.inFile
    EEG = pop_loadset(inFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select only fixations to faces
if filter_faces
    event_faces = [];
    for ev = 1:numel([EEG.event]) 
        %keyboard
        if EEG.event(ev).faces == 1
            event_faces = [event_faces ev];
        end
    end
    EEG.event(event_faces) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select only fixations in the exploration trials
if filter_EX
    event_EX = [];
    for ev = 1:numel([EEG.event]) 
        %keyboard
        if EEG.event(ev).VSNT_EX == 0
            event_EX = [event_EX ev];
        end
    end
    EEG.event(event_EX) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select only fixations in the active search trials
if filter_AS
    event_AS = [];
    for ev = 1:numel([EEG.event]) 
        %keyboard
        if EEG.event(ev).VSNT_EX == 1
            event_AS = [event_AS ev];
        end
    end
    EEG.event(event_AS) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tend_load = toc(tstart_load)
    
%     latency = [EEG.event.latency];
%     perm = randperm(length(latency));
%     PermLatency = latency(perm);
%     for il = 1:length(EEG.event)
%         EEG.event(il).OrigLatency = latency(il)
%         EEG.event(il).PermLatency = PermLatency(il)
%     end

    %DEFINING THE DESIGN
    
    % tstart_unfold = tic;
    
    cfgDesign = [];
    cfgDesign.eventtypes = {'fixation'}; % we model the fixation onset
    %cfgDesign.formula = {'y ~ 1+cat(isEasy)'} %cat(stimulusType)'}%+cat(fixationRank)'};
    cfgDesign.formula = model ; %saccade_amp + cat(hard_easy)'}%+stimulusDur+cat(hard_easy)'}
    
    %One needs to be careful to not overfit the data, regularization or cross validation can help here.
    EEG = uf_designmat(EEG,cfgDesign);
    
    % DELETE NANS mdf 29/6/21 
    cfgNueva = [];
    cfgNueva.method = 'drop';
    EEG = uf_imputeMissing(EEG,cfgNueva)
    
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

%     for il = 1:length(EEG.event)
%         EEG.event(il).latency = EEG.event(il).PermLatency 
%     end   
    EEG.event(ismember({EEG.event.type},{'bad_ET','boundary'})) = []; %DAC mod 30 sept 2021

    %FITTING MODEL 
    EEG= uf_glmfit(EEG);
    
    % tend_unfold = toc(tstart_unfold)
    
    % tstart_save = tic;

    suje = names{iSuj}; sujName = suje(1:3);
    
    n = struct(cond{1},sum([EEG.event.(hipothesis)]==1),cond{2},sum([EEG.event.(hipothesis)]==0));
    EEG.unfold.N = struct(cond{1}, n.(cond{1}), cond{2},n.(cond{2})); % unreliable
    
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
    
    % tend_save = toc(tstart_save)
    
end

%% 3rd step Generate Data matrix for TFCE and plots for variables of interest

% unfiltered data has not the baseline correction

inFolder                =  fullfile(session_path.out, 'Unfold2021/' );
models = fullfile(inFolder, 'models/' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use the name of the folder where the model were stored
model2analyse = modelName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if freqFilteredEEG
    folderOut = fullfile(models,currentBand.name, '/',model2analyse);
else   
    folderOut = fullfile(models,model2analyse);
end

if exist(folderOut, 'dir')
    
    load(fullfile(folderOut, 'Results.mat'));
    variableStruct =[];
    subjects   = fieldnames(unfoldResults);
    for iSuj  = 1:length(fieldnames(unfoldResults))
        suj   = subjects{iSuj}; 
        result = unfoldResults.(suj);
        variableStruct  = fun_makeStructToERPunfold(variableStruct,variablesOI,iSuj,suj,result,hipothesis)
    end

    fileOut   = fullfile(folderOut,'variableStruct.mat');

    save(fileOut,'variableStruct')

%load(fileOut)
else
    fprintf("the model doesn't  exist, check model name or fit it \n")
end

%% 4th step run TFCE
variable2TFCE   = 'NCRank';%'(Intercept)';

inFolder                =  fullfile(session_path.out, 'Unfold2021/' );
models = fullfile(inFolder, 'models/' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use the name of the folder where the model were stored
model2analyse = modelName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if freqFilteredEEG
    folderOut = fullfile(models,currentBand.name, '/',model2analyse);
else   
    folderOut = fullfile(models,model2analyse);
end

fileStruct  = fullfile(folderOut,'variableStruct.mat');
if exist(folderOut, 'dir')
    load(fileStruct);
    times   = variableStruct.times; 
    chanlocs = load(fullfile(session_path.aux_data,'locationFile.mat'));
    chanlocs = chanlocs.e_loc;
    if baselineCorr == 1 && ~freqFilteredEEG
        variableStruct = fun_applyBaselineCorrection(variableStruct,variable2TFCE,freqFilteredEEG)
        fileOut   = fullfile(folderOut,'beta_', variable2TFCE,'_BC.mat');

        save(fileOut,'variableStruct')    
    end
    tfceResult  = fun_runTFCE(fp,folderOut,variableStruct, variable2TFCE,chanlocs)
else
    fprintf("the model doesn't  exist, check model name and variableStruct \n")
end
 
%% Functions 

function ERPstruct  = fun_makeStructToERPunfold(ERPstruct,variablesOI,iSuj,suj,result,hipothesis)
    

    ERPstruct(iSuj).times    = result.times;
    ERPstruct(iSuj).Suj      = suj;
%     if sum(contains(variablesOI, hipothesis))>0
%         cond = split(hipothesis,'_');
%         ERPstruct(iSuj).(['n_' cond{1}]) = result.N.(cond{1});
%         ERPstruct(iSuj).(['n_' cond{2}]) = result.N.(cond{2});
%     end

    for var = 1:length(variablesOI)
        covariate             = find(ismember(result.variablenames,variablesOI{var}));
    %                     if sum(covariate)~=1
    %                         fprintf('unmatch or repeted variable')
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
