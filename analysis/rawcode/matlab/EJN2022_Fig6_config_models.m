%Edit this lines to choose a model to fit It should be written in Wilkinson
%notation, categorical and variables model with splines should be indicated
%cat(categorical variable), spl(variable,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fig. 4
% % fixs_csv_file  = '16subj_VSNTpostvsEX_for_unfold_GlobalStandarisedCenteredRank.csv';% OLD This is the file used in EJN Round 0 for VSNTpre_EX condition. It has 5127 fixs.
% % fixs_csv_file  = 'EJN2022_16subj_VSNT_EX.csv';                                      % NEW This is the file used in EJN Round 1 for VSNT_EX condition. It has 6421 fixs.
% % fixs_csv_file  = 'EJN2022_16subj_VSNTpre_EX.csv';                                   % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% fixs_csv_file  = 'EJN2022_16subj_VSNTpre_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'pre_';
% 
% hipothesis      = 'faces';
% modelName       = 'Inter+face';
% model           = {'y ~ 1 + cat(faces)'}; %+ cat(VSNT_EX) + StandFixRank
% 
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% variablesOI     = {'(Intercept)' hipothesis};
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% condi            = {'O','F'};%%%condi = split(hipothesis,'_');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fig. 5
% fixs_csv_file  = 'EJN2022_16subj_VSNTpre_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'pre_';
% 
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+face';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces)'}; %+ cat(VSNT_EX) + StandFixRank
% 
% variablesOI     = {'(Intercept)' hipothesis 'faces'};
% condi            = {'EX','AS'};%%%condi = split(hipothesis,'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fig. 6
% fixs_csv_file  = 'EJN2022_16subj_VSNTpre_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'pre_';
% 
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+faces+NCrank+interTaskNCRank';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank
% 
% variablesOI     = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'};
% condi            = {'EX','AS'};%%%condi = split(hipothesis,'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % New analysis: Interaction between Task and Category (same data)
% fixs_csv_file   = 'EJN2022_16subj_VSNTpre_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'pre_';
% 
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+faces+interTaskFaces';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + faces:VSNT_EX'};
% 
% variablesOI     = {'(Intercept)' hipothesis 'faces' 'VSNT_EX:faces'};
% condi            = {'EX','AS'};%%%condi = split(hipothesis,'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % New analysis: Full model like figure 6 but with absent trials instead pretarget
% fixs_csv_file   = 'EJN2022_16subj_VSNTabsent_EX_manyNranks.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'abs_';
% 
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+faces+NCrank+interTaskNCRank';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank
% 
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% variablesOI     = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'};
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% condi            = {'EX','AS'};%%%condi = split(hipothesis,'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % New analysis: Full model like figure 6 but with absent trials instead
% % exploration
% 
% fixs_csv_file   = 'EJN2022_16subj_VSNTpretarget_absent.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'abspre_';
% 
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+faces+NCrank+interTaskNCRank';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank
% 
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% variablesOI     = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'};
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% condi            = {'EX','AS'};%%%condi = split(hipothesis,'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % New analysis: Full model like figure 6 but with (short) absent trials instead pretarget
% fixs_csv_file   = 'EJN2022_16subj_VSNTabsentshort_EX.csv';                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
% datasetID       = 'absshort_';
% 
% hipothesis      = 'VSNT_EX';
% modelName       = 'Inter+task+faces+NCrank+interTaskNCRank';
% model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank
% 
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% variablesOI     = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'};
% %e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
% condi            = {'EX','AS'};%%%condi = split(hipothesis,'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New analysis: Full model like figure 6 but with (short) absent trials instead pretarget
if ~exist('RR')
    fprintf('ERROR: You are running the Absent step analysis and RR (max rank allowed) is not defined.\n')
end
fprintf('******************************************\n')
fprintf('\t\tMoving rank: %d\n',RR)
fprintf('******************************************\n')
fixs_csv_file   = sprintf('EJN2022_16subj_VSNTabsentshort_EX_%d.csv',RR);                          % NEW This is the file used in EJN Round 1 for VSNTpre_EX condition. It has 5578 fixs.
datasetID       = sprintf('absshort_%d_',RR);

hipothesis      = 'VSNT_EX';
modelName       = 'Inter+task+faces+NCrank+interTaskNCRank';
model           = {'y ~ 1 + cat(VSNT_EX) + cat(faces) + NCRank + NCRank:VSNT_EX'}; %+ cat(VSNT_EX) + StandFixRank

%e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
variablesOI     = {'(Intercept)' hipothesis 'faces' 'NCRank' 'NCRank:VSNT_EX'};
%e.g.: '(Intercept)' , 'hard_easy', 'stimulusDur', 'saccade_amp'
condi            = {'EX','AS'};%%%condi = split(hipothesis,'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('CSV file: %s\n', fixs_csv_file);
fprintf('Models name: %s\n', modelName);
fprintf('Model (Wilkinson notation): %s\n', model{1});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baselineCorr    = 1;
% ourPermTest     = 1;
fprintf('Analysing with baseline correction\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter conditions % jek 10/11/2022
filter_faces    = 0; % 1 if filter fix to faces O otherwise
filter_EX       = 0;
filter_AS       = 0;