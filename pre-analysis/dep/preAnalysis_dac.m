%% Settings
clear all
close all
%cd /Users/dac/Documents/GitLab/visualsearch_eegem/codes

cd /Users/dac/Documents/GitLab/corregistro/Analysis_dac%/media/dac/data/repos/corregistro/codes
[runpath,code_path,session_path,whoisrunning]=add_paths_EEG_ET_beta(1,1)%(1,1) for use in laptop

eeglabpath = code_path.eeglabpath
%desktop agregar a addppath func
%addpath('~/Descargas/codes')
addpath(eeglabpath); 

ftpath      =        code_path.toolbox ;

addpath(ftpath); 
addpath(code_path.corregistro_fn)
addpath(code_path.my_functions)
clear fp; fp = FP_basePreAnalysis('dac');

eeglab

datapath =  code_path.raw;

cd(datapath)
%% STEP 1 EEG preprocessing
if 1
folderIn =   datapath;
folderOut =  session_path.renamed;
namesIn  = fp.fn.loadSubjects(folderIn, 1, 1,'bdf','EEGEYE');
namesOut = fp.fn.loadSubjects(folderOut, 1, 1, 'set');


for iSuj = 1:length(namesIn)
   suj = namesIn{iSuj}
    if ~ismember([suj '.set'], namesOut)
        file = fullfile(folderIn, suj);
        % Load
        fp.cfg.inFile  = [ folderIn   suj ];
        fp.cfg.outFile = [ folderOut  suj '.set']; 
        fp.cfg.ref     = [65 66];  % vector of references for pop_biosig
        fp.cfg.nChansKeep = 64;     % number of chans with data 64

        EEG = fp.EEG_preprocessing();
        %AVERAGE
        EEG = pop_reref( EEG, []);
        EEG = pop_saveset(EEG, fp.cfg.outFile);
    end
end
end
%% STEP 2  EEG preFiltering
folderIn  = session_path.renamed;
folderOut = session_path.prefiltered;
namesIn  = fp.fn.loadSubjects(folderIn, 1, 1, 'set');
namesOut = fp.fn.loadSubjects(folderOut, 1, 1, 'set');
fp.cfg.chanlocsFilePath = [code_path.eeglabpath '/functions/resources/Standard-10-5-Cap385.sfp' ]; %complete path to chanloc file

for iSuj = 1:length(namesIn)
    suj = namesIn{iSuj}; suj = suj(1:end-4);
    if ~ismember([suj,'.set'], namesOut)
        fp.cfg.inFile  = [ folderIn   suj '.set'];
        fp.cfg.outFile = [ folderOut  suj '.set']; 
        fp.cfg.preFilterEdges   = [.1 100]; % edges for band pass filter using pop_eegfiltnew()
        fp.cfg.noFiltAnalog     = 0; % [Bool] wheter to remove analog before filtering (1) or not (0)
        
        EEG = fp.EEG_prefiltering();
        EEG = pop_saveset(EEG, fp.cfg.outFile);
    end
end
%% STEP 3.1 Slect channels to interpolate
folderIn = session_path.prefiltered;
names    = fp.fn.loadSubjects(folderIn, 1, 1, 'set'); 

for iSuj = 1% 1:length(names)
    suj = names{iSuj};
    fp.cfg.inFile = [folderIn suj];
    fp.cfg.nHeadChans = 64;
    EEG = fp.EEG_selectChanInterpolate()
    keyboard
    close all
end
%% STEP 3.2 Interpolate
% Choose between 1 to 64

% Notes: recheck 2020 channels inspected are marked with **
                
                %21 noviembre 2019----
badchanslist = {'EEGEYEE01', [52];%<--average[]; ... visual:  , espectro: 52 , variance:ok----** mastoide?       
                'EEGEYEE02', [28 2];%<--average[27,28]; ... visual:ok, espectro: 27,28, variance:ok---- **
                'EEGEYEE03', [34 52 54];%<--average[19, 20, 26]; ... visual: {19. 20. 21,30,31,32,56,57,58}, espectro y var : 19 , 20 , 26 **                                            
                'EEGEYEE04', [24];%<--average...[3, 35, 42];...  visual: {} small amplitude and noisy  (periods with bad signals) , espectro: {2}(pico muy poco definido a 10hz),variance:**malo
                %21 noviembre 2019
                'EEGEYEE05', [ 34, 35];...      visual:{34,35} , espectro: {34, 35,41},variance: ok---- **
                'EEGEYEE06', [34, 35];...        visual: {34,35}, espectro: 35,variance: ok---- ** 
                'EEGEYEE07', [33 34 35];...  visual: {}, espectro: 33 34 35, variance:ok----    **
                'EEGEYEE08', [];...          visual: 1 Very nosy subject , secments with high amplitude variantion, espectro: shows low amplitud for 10hz region,variance:ok **
                'EEGEYEE09', [];...      visual: {28} (O2 POz Iz O1 seems to have some high freq noise, not sure) , espectro: {28}it has a peak near 26hz for several channels,variance:ok **
                'EEGEYEE10', [46];%<--average;[24 ];...      visual:, espectro: {46} small response for all channels at 10hz specially cannel 46,variance:35 bit high ** 
                'EEGEYEE11', [34 35 42];%average...[];...        visual: {1}, espectro: {7,9,35,42} abnormal spectral distribution,variance: 7 ,42 high(over 10⁴) ** mastoide?
                'EEGEYEE12', [52 57];%<--average;[52];...        visual: {52}, espectro: {2,52},variance: 1,2,33,34,35 **
                'EEGEYEE13', [24 42 ];... visual: 42, espectro: 42, variance: ok **
                'EEGEYEE14', [24];%<--average;[];... visual: , espectro: 16,23,24,28,variance: 1,35---- **
                'JK_EEGEYE4', [55];... visual: , espectro: 55,variance: 40,41,55,62----  nbch80 ** canales con ruido como que se desconectan y poca duracion
                'Matias_EEGEYE', [28,34];... visual: , espectro: 28,34,variance: 34----     ** eventos solo del tipo 49 !!!!!!
                'EEGEYEJ01', [2 35];... visual: , espectro: ,variance:      nbch-- correr denuevo correr denuevo correr denuevo ** mastoide?-con average tampoco mejora
                'EEGEYEJ02', [21,63];... visual: 21 , espectro:21 ,63 ,variance:21,63     ** not enough events just 23 !!!!!!
                'EEGEYEJ03', [42];... visual: 35, espectro: 42 ,variance: 42   **
                'EEGEYEJ04', [];... visual: , espectro: 14,9,44,51 ,variance:    **
                'EEGEYEJ05', [35 63];... visual: , espectro: 35,63 ,variance: 35 ,63 **
                'EEGEYEJ06', [];... visual: 15,52 , espectro: 61,15,52,variance: ok  **
                'EEGEYEJ07', [15];... visual:15 , espectro:15,26,variance: 26   **
                'EEGEYEJ08', [];... visual:24,25,28,64 , espectro:,variance:   nbch-- Peak at appr.26hz for all channels **
                'EEGEYEJ09', [63];... visual:34,35,52 , espectro:,variance:26,57,63    nbch-- high noise in several channels from 15 hz and beyond.**
                'EEGEYEJ10', [20,35,63];... visual: 63 , espectro:20,35,variance:20  nbch-- very small 10hz amplitud **
                'EEGEYEJ11', [];... visual: 8,9 , espectro:6,9,44 variance:ok  **
                'EEGEYEJ12', [];... visual: , espectro: variance:  nbch-- about 15hz peak **
                
             
                };
folderIn            =  session_path.prefiltered;
folderOut           =  session_path.interpolated ;
names               =  fp.fn.loadSubjects(folderIn, 1, 1, 'set');
fp.cfg.badchanslist =  badchanslist;
for iSuj=15%:length(names)
    suj = names{iSuj};
    fp.cfg.inFile = [folderIn suj];    
    EEG = fp.EEG_interpolate();
    EEG = pop_saveset(EEG,[ folderOut  suj]);          
end




%% STEP 4.1 ET Preprocessing
folder    =  session_path.interpolated ;% EEG files to be synchronized
folderOut = session_path.et_mat;
names     =  fp.fn.loadSubjects(folder, 1, 1, 'set');
namesOut = fp.fn.loadSubjects(folderOut, 1, 1, 'mat');
FolderIn  =   session_path.rawet;
% RESOLVER ERROR EN ET DE J02 
for iSuj = 1:length(names)
    
    filename = names(iSuj);
    suj      = filename{1}(1:end-4);
    if ~ismember([suj,'_ET_EYEMAP_events.mat'], namesOut)
        ind      = find(strcmp(session_path.subjname, suj(end-2:end)));

        if ~isempty(ind)
            fp.cfg.inFile = [FolderIn '/'  session_path.sessionfilenames{ind} '.asc'];
            fp.cfg.outFile = [session_path.et_mat  suj '_ET_EYEMAP_events.mat']
            ET =  fp.ET_preprocessing()
            %keyboard
        end
    end
end



%% STEP 4.2 Sync EEG and ET   
%[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_EEG(1,[]);%leave usr=[]

folderIn    = session_path.interpolated;
folderOut   = session_path.synced;  
etMATfolder = session_path.et_mat; 
names       = fp.fn.loadSubjects(folderIn, 1, 1, 'set');
namesOut    = fp.fn.loadSubjects(folderOut, 1, 1, 'set');



%sujeto E3 falta el ultimo eyemap 
%J1(13) solo la prim eymap y J2(14) no
%sincronizan J04 no sincroniza todo los puntos J5 no sincroniza eyemap del
%medio J7 J8 J9 no sincroniza en el medio 


for iSuj = [17 20]%:length(names)
    close all
    filename = names(iSuj);
    suj      = filename{1}(1:end-4);
    ind      = find(strcmp(session_path.subjname, suj(end-2:end)));
    if 1%~ismember([suj '_sync.set'], namesOut)
        if ~isempty(ind)
            fp.cfg.inFileEEG = [folderIn suj '.set'];
            fp.cfg.inFileET  = [ etMATfolder suj '_ET_EYEMAP_events.mat'];
            EEG              = fp.EEG_ET_sync()
            keyboard
            EEG = pop_saveset(EEG,[ folderOut  suj '_sync.set']);
        end
    end
% eeglab
end
%% STEP 5 Engbert

folderIn    = session_path.synced;
folderOut   = session_path.eng;  
names       = fp.fn.loadSubjects(folderIn, 1, 1, 'set');

%!!!E06esta raro

for iSuj =6%1:length(names)
    close all
    filename      = names(iSuj);
    suj           = filename{1}(1:end-4);
    fp.cfg.inFile = [folderIn suj '.set'];

    if 1%~isempty(ind)
        try
            fp.cfg.eye     = 'R';%eye to be used 'R' 
            EEG           = fp.engbertDetection()        
        catch
            fprintf('changing eye selection to L \n')
            fp.cfg.eye     = 'L';%eye to be used  'L'
            EEG           = fp.engbertDetection() 
        end
%         fp.cfg.eye                  = 'R';%eye to be used 'R' or 'L'
%         fp.cfg.inFile = [folderIn suj '.set'];
%         EEG              = fp.engbertDetection()
        keyboard
        EEG = pop_saveset(EEG,[ folderOut  suj '_eng.set']);
    end
    
end
%% (incomplete)!!!!!STEP 5.1 EOG Surrogates
%I'm doing this after having tested the whole pipeline without step 5.1 ,
%since I don't have EOG channels in my data I want to see if
%performance of the artifact rejection algorithm can be improved adding
%rescale eyetracker data as fake EOG channels


folderIn    = session_path.eng;
folderOut   = session_path.eog;  
names       = fp.fn.loadSubjects(folderIn, 1, 1, 'set');



for iSuj =1:length(names)
    close all
    filename      = names(iSuj);
    suj           = filename{1}(1:end-4);
    fp.cfg.inFile = [folderIn suj '.set'];

    if 1%~isempty(ind)
        EEG                        = pop_loadset(inFile); 

        EEG = pop_saveset(EEG,[ folderOut  suj '.set']);
    end
    
end


%% STEP 6 ICA Training

folderIn    = session_path.eng;
folderOut   = session_path.icaWeights;  
names       = fp.fn.loadSubjects(folderIn, 1, 1, 'set');
namesOut       = fp.fn.loadSubjects(folderOut, 1, 1, 'mat');


for iSuj = 21%1:length(names)
 
    
    close all
    filename      = names(iSuj);
    suj         = filename{1}(1:end-12);
    fp.cfg.inFile = fullfile(folderIn, filename);
    fp.cfg.sacLabel = 'saccade';
    fp.cfg.fixLabel = 'fixation';

    
   
    
    
    if ~ismember([suj 'weights.mat'], namesOut)
            [sph wts]   = fp.trainICA()
            weights.sph = sph;
            weights.wts = wts;
            %get subject name and save
            fprintf('saving struct with weights for subject %s', suj)
            %output file must be complete path and .mat file
            fileOut    = [folderOut suj 'weights.mat'];
            save(fileOut,'-struct','weights')


        
    end
end


%% STEP 7 epoching data


folderIn    =  session_path.eng;
folderOut   =  session_path.epoched;  
names       =  fp.fn.loadSubjects(folderIn, 1, 1, 'set');
namesOut       =  fp.fn.loadSubjects(folderOut, 1, 1, 'set');



for iSuj = 21%1:length(names)
    close all
    filename      = names(iSuj);
    fp.cfg.inFile = fullfile(folderIn ,filename);
    %fp.cfg.chanlocsFilePath = '/Users/dac/Documents/MATLAB/eeglab14_1_2b/functions/resources/Standard-10-5-Cap385.sfp' ;
    
    
    
    suj         = filename{1}(1:end-12);

    if ~ismember( [suj 'stim_epoched.set'], namesOut)
        EEG   = fp.epochData()
       
        %get subject name and save
        fprintf('saving epoched EEG for subject %s', suj)
        %output file must be complete path and .mat file
        fileOut    = [folderOut  suj 'stim_epoched'];
        EEG = pop_saveset(EEG,[ fileOut '.set']);
    
        
    end
    
end



%% STEP 8 IC's selection

folderIn      =  session_path.eng;
folderEpoch   =  session_path.epoched;
folderWeights =  session_path.icaWeights; 
folderOut     =  session_path.data_analysis;  
names         =  fp.fn.loadSubjects(folderIn, 1, 1, 'set');
namesEpoch    =  fp.fn.loadSubjects(folderEpoch, 1, 1, 'set');%make check for non existent files later
namesOut       =  fp.fn.loadSubjects(folderOut, 1, 1, 'set');

fp.cfg.sacstring = 'saccade';
fp.cfg.fixstring = 'fixation';
fp.cfg.plot      = 1;
for iSuj =21 %1:length(names)
    close all
    filename            = names(iSuj);
    fp.cfg.inFile       = fullfile(folderIn, filename);
    suj                 = filename{1}(1:end-12);
    %ind                 = find(ismember(names,namesEpoch))  

    if  1%~ismember( [suj 'analysis.set'], namesOut)
        fp.cfg.weights      = [folderWeights suj 'weights.mat'];
       
        EEG     = fp.selectIcPlochl()
        keyboard
        fp.cfg.badcomp      = EEG.reject.gcompreject;
        
        fprintf('Loading epoched EEG for subject %s for removing artifact related channels\n', suj)
        epocheFileIn        = [folderEpoch  suj 'stim_epoched.set'];
        fp.cfg.inFile       = epocheFileIn; % epoched data file to be plöcheled
        EEG                 = fp.removeComp()
        fileOut    = [folderOut suj 'analysis.set'];
        EEG = pop_saveset(EEG,fileOut );
    
        
    end
    
end


