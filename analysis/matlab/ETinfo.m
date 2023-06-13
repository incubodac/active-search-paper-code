%% mdf 18/3/21 

clear all
close all
cd /home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/functions
%cd /media/dac/data/repos/corregistro/codes
[runpath,code_path,session_path,whoisrunning]=add_paths_EEG_ET_beta(1,2);%(1,1) for use in laptop

eeglabpath = code_path.eeglabpath;
addpath(eeglabpath); 

ftpath      =        code_path.toolbox ;
%addpath('/media/dac/data/repos/visualsearch_eegem/codes')
addpath(ftpath); 
addpath(genpath(code_path.corregistro_fn))
addpath(code_path.my_functions)
addpath(code_path.psyco )
addpath(runpath);
addpath(code_path.functions)
clear fp; fp = FP_epochAnalysis();

eeglab
matdir              = session_path.matfiles;
datapath            = session_path.data_analysis; 
cd(datapath)

inFolder                =  session_path.data_analysis;
fp.cfg.subjname         =  session_path.subjname;
fp.cfg.sessionfilenames =  session_path.sessionfilenames;
fp.cfg.matdir           =  matdir;
folderOut               =  session_path.conditions_data ;

names = fp.loadSubjects(inFolder, 1, 1, 'set');

%% export fixs struct for every subject

inFolder                =  session_path.data_analysis;
fp.cfg.subjname         =  session_path.subjname;
fp.cfg.sessionfilenames = session_path.sessionfilenames;
fp.cfg.matdir           = matdir;
% folderOut               = session_path.fixs ;
% TOTEST
folderOut               = '/media/cbclab/MARIADAFON1T/Analysis_2020/fixs2/';

names = fp.loadSubjects(inFolder, 1, 1, 'set');
minDur4fix = 0.1;
for iSuj = [1:4 6:11 13:15 19 22:23]%subjects with correct scanpath%[1:4 6:11 15:17]%[12:15]%15:length(names)

    % mdf 26/4/21 avoid iSuj 21 which is J10 participant
    
    fp.cfg.inFile = [inFolder names{iSuj}];

    [EEG fixs trials ExpTrials]   = fp.fixationsAnalysisOne(minDur4fix,0,1) %Here is when mat files for each subject  is required Joe series problem.

    % [EEG fixs trials ExpTrials]   = fp.fixationsAnalysisAll() %Here is when mat files for each subject  is required Joe series problem.

    suj = names{iSuj}; sujName = suj(7:end-13);

    fp.cfg.outFile = [ folderOut  sujName '_fixEpoch.set' ]; 
    
    fileOut = [folderOut  sujName 'FIXstruct.csv'];

    writetable(struct2table(fixs), fileOut)
%26/02 epocheado para distractores en todas las condiciones caras vs
%objetos y para targets solo VS correctos  detalles en (indexCondEpoch)

end