%% mdf 5/3/21 

clear all
close all
cd /home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/functions
%cd /media/dac/data/repos/corregistro/codes
[runpath,code_path,session_path,whoisrunning]=add_paths_EEG_ET_beta(1,2);%(1,1) for use in laptop

eeglabpath = code_path.eeglabpath;
addpath(genpath(eeglabpath),'-end');

ftpath      =        code_path.toolbox ;
%addpath('/media/dac/data/repos/visualsearch_eegem/codes')
addpath(ftpath); 
addpath(genpath(code_path.corregistro_fn))
addpath(code_path.my_functions)
addpath(code_path.psyco )
addpath(runpath);
addpath(code_path.functions)
clear fp; fp = FP_epochAnalysis();

matdir              = session_path.matfiles;
datapath            = session_path.data_analysis; %session_path.data_analysis500;
cd(datapath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inFolder                =  session_path.data_analysis; %session_path.data_analysis500; 
% inFolder                =  session_path.BrunoData;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp.cfg.subjname         =  session_path.subjname;
fp.cfg.sessionfilenames =  session_path.sessionfilenames;
fp.cfg.matdir           =  matdir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderOut               =  session_path.freqEEG ;
% folderOut               =  fullfile(session_path.out,'/brunobian/freqEEG/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = fp.loadSubjects(inFolder, 1, 1, 'set');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subj_index = [1:4 6:11 13:15 19 22:23]; % mdf 13/5/21 sobre qué sujetos epochear, etc
% subj_index = [1:length(names)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose the freq-band to filter
Theta = struct('limits',[4  8],'name','Theta');
Alpha = struct('limits',[8 13],'name','Alpha');
Beta = struct('limits',[15 25],'name','Beta');
LowBeta = struct('limits',[14 19],'name','LowBeta');
HighBeta = struct('limits',[20 30],'name','HighBeta');
LowGamma  = struct('limits',[30 40],'name','LowGamma');

currentBand = LowGamma;

eeglab

%% resulting EEG after filtering and applied the hilbert transform 

minDur4fix = 0.1;
for iSuj = subj_index

    fp.cfg.inFile = [inFolder names{iSuj}];

% filter by freq band
    EEG = pop_loadset(fp.cfg.inFile);
    EEG = eeg_checkset(EEG);
    EEG = pop_eegfiltnew(EEG, 'locutoff',currentBand.limits(1),'hicutoff', currentBand.limits(2),'minphase',true);
    for ch=1:size(EEG.data,1)
        EEG.data(ch,:,:) = imag(hilbert(EEG.data(ch,:,:)));
    end
    
% filtering in special freqs    
%     for min_freq = 3:30
%         EEG   = fp.filterband(fp.cfg.inFile, min_freq, min_freq+1,true);
% 
        suj = names{iSuj}; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sujName = suj(7:end-13); % filter by band
        % sujName = suj(1:end-4); %filter by special freq
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fp.cfg.outFile = [ folderOut  currentBand.name '/' sujName '_freq_' currentBand.name '.set' ]; % save filter by band
        %fp.cfg.outFile = [ folderOut  'multipleFreq500/' sujName '_freq_'
        %num2str(min_freq) '.set' ]; % save filter by special freq
        % fp.cfg.outFile = [ '/media/cbclab/MARIADAFON1T/Analysis_2020/'  sujName '_new_' currentBand.name '.set' ]; 

        pop_saveset(EEG, fp.cfg.outFile);
end


%% freq vs time
folderIn = [folderOut 'multipleFreq500/'];

names = {'E01', 'E03', 'E04', 'E05', 'E07','E08', 'E09', 'E10',  'E11', 'E12', 'E14', 'J03', 'J04', 'J08', 'J11', 'J12'};

type_of_trial = readtable([folderIn 'VSac_EX.csv']);
VSac_trials = type_of_trial(type_of_trial.VSac_EX==1,:); % a=absent; c=correct
EX_trials = type_of_trial(type_of_trial.VSac_EX==0,:);

type_of_trial2 = readtable([folderIn 'VSpcs_EX.csv']);
VSpcs_trials = type_of_trial2(type_of_trial2.VSpcs_EX==1,:); % =present; c=correct; s=seen (at least one fix to target longer than 200 ms)

freqStruct = []; LinRegrTime = [];
VSac_Tsubj = []; EX_Tsubj = []; VSpcs_Tsubj = []; VSpcsNan_Tsubj = []; VSac_menos_EX_Tsubj = []; 
VSpcs_menos_EX_Tsubj = []; VSpcs_menos_VSac_Tsubj = []; 
for iSuj = 1:2%length(names)%subj_index 
    VSac_subj = VSac_trials(strcmp(VSac_trials.trial_subject, char(names(iSuj))),:);
    EX_subj = EX_trials(strcmp(EX_trials.trial_subject, char(names(iSuj))),:);
    VSpcs_subj = VSpcs_trials(strcmp(VSpcs_trials.trial_subject, char(names(iSuj))),:);
    
    if strcmp('J08', char(names(iSuj)))
        VSac_subj = VSac_subj(VSac_subj.trial_number<166,:);
    end
    VSac_freq = []; EX_freq = []; VSacmenosEX_freq = []; VSpcs_freq = []; VSpcsmenosEX_freq = []; VSpcsmenosVSac_freq = [];
    VSpcsNan_freq = [];
    
    B1_VSpcsNan_freq = []; B2_VSpcsNan_freq = []; R2_VSpcsNan_freq = []; p_VSpcsNan_freq = [];
    for min_freq = 3:5%30
            % inFile = [folderIn 'EEGEYE' char(names(iSuj)) '_analysis_freq_' num2str(min_freq) '.set'];
            inFile = [folderIn 'EEGEYE' char(names(iSuj)) '_analysis-500_freq_' num2str(min_freq) '.set'];
            
            EEG = pop_loadset(inFile);
            if min_freq == 3 && iSuj==1
                freqStruct.times = EEG.times;
                freqStruct.chanlocs = EEG.chanlocs(1:64);
                channels = {EEG.chanlocs.labels};
                times = EEG.times;
            end
            
            X = [ones(size(times')) times'];
            
            data = EEG.data;
            data = fun_baselineCorrection(data,times);

            VSac_freq = fun_freqChanTime(data,VSac_subj,VSac_freq);

            EX_freq = fun_freqChanTime(data,EX_subj,EX_freq);
            
            VSpcs_freq = fun_freqChanTime(data,VSpcs_subj,VSpcs_freq);
            
            VSpcsNan_freq = fun_freqChanTime_nan(data,VSpcs_subj,VSpcs_freq,times);
            
            % linear regressions in time
            B1_VSpcsNan_ch = []; B2_VSpcsNan_ch = []; R2_VSpcsNan_ch = []; p_VSpcsNan_ch = [];
            for ch = 1:length(channels)   
                [b, bint, ~, ~, stats] = regress(VSpcsNan_freq(ch,:,min_freq-2)',X);
                B1_VSpcsNan_ch = [B1_VSpcsNan_ch b(1)];
                B2_VSpcsNan_ch = [B2_VSpcsNan_ch b(2)];
                R2_VSpcsNan_ch = [R2_VSpcsNan_ch stats(1)];
                p_VSpcsNan_ch = [p_VSpcsNan_ch stats(3)];
            end
            B1_VSpcsNan_freq = [B1_VSpcsNan_freq; B1_VSpcsNan_ch];
            B2_VSpcsNan_freq = [B2_VSpcsNan_freq; B2_VSpcsNan_ch];
            R2_VSpcsNan_freq = [R2_VSpcsNan_freq; R2_VSpcsNan_ch];
            p_VSpcsNan_freq = [p_VSpcsNan_freq; p_VSpcsNan_ch];
    end
    LinRegrTime.VSpcsNan.(char(names(iSuj))).B1 = B1_VSpcsNan_freq;
    LinRegrTime.VSpcsNan.(char(names(iSuj))).B2 = B2_VSpcsNan_freq;
    LinRegrTime.VSpcsNan.(char(names(iSuj))).R2 = R2_VSpcsNan_freq;
    LinRegrTime.VSpcsNan.(char(names(iSuj))).p = p_VSpcsNan_freq;
    
    VSac_Tsubj = cat(4,VSac_Tsubj,VSac_freq);
    EX_Tsubj = cat(4,EX_Tsubj,EX_freq);  
    VSpcs_Tsubj = cat(4,VSpcs_Tsubj,VSpcs_freq); 
    VSpcsNan_Tsubj = cat(4,VSpcsNan_Tsubj,VSpcsNan_freq); 
    VSac_menos_EX_Tsubj = cat(4,VSac_menos_EX_Tsubj,VSac_freq-EX_freq);
    VSpcs_menos_EX_Tsubj = cat(4,VSpcs_menos_EX_Tsubj,VSpcs_freq-EX_freq);   
    VSpcs_menos_VSac_Tsubj = cat(4,VSpcs_menos_VSac_Tsubj,VSpcs_freq-VSac_freq); 
%     
%     freqStruct(iSuj).('VSac') = VSac_freq;    
%     freqStruct(iSuj).('EX') = EX_freq;  
%     freqStruct(iSuj).('VSpcs') = VSpcs_freq;     
%     freqStruct(iSuj).('VSpcsNan') = VSpcsNan_freq;     
%     freqStruct(iSuj).('VSac_menos_EX') = VSac_freq-EX_freq;   
%     freqStruct(iSuj).('VSpcs_menos_EX') = VSpcs_freq-EX_freq;       
%     freqStruct(iSuj).('VSpcs_menos_VSac') = VSpcs_freq-VSac_freq; 
end

freqStruct.VSac.mean = nanmean(VSac_Tsubj,4);
% freqStruct.VSac.std = std(VSac_Tsubj(~isnan(VSac_Tsubj)),4);
    
freqStruct.EX.mean = nanmean(EX_Tsubj,4);  
% freqStruct.EX.std = std(EX_Tsubj(~isnan(EX_Tsubj)),4);  

freqStruct.VSpcs.mean = nanmean(VSpcs_Tsubj,4); 
% freqStruct.VSpcs.std = std(VSpcs_Tsubj(~isnan(VSpcs_Tsubj)),4); 
    
freqStruct.VSpcsNan.mean = nanmean(VSpcsNan_Tsubj,4); 
% freqStruct.VSpcsNan.std = std(VSpcsNan_Tsubj(~isnan(VSpcsNan_Tsubj)),4); 
    
freqStruct.VSac_menos_EX.mean = nanmean(VSac_menos_EX_Tsubj,4);
% freqStruct.VSac_menos_EX.std = std(VSac_menos_EX_Tsubj(~isnan(VSac_menos_EX_Tsubj)),4);
    
freqStruct.VSpcs_menos_EX.mean = nanmean(VSpcs_menos_EX_Tsubj,4);   
% freqStruct.VSpcs_menos_EX.std = std(VSpcs_menos_EX_Tsubj(~isnan(VSpcs_menos_EX_Tsubj)),4);   
    
freqStruct.VSpcs_menos_VSac.mean = nanmean(VSpcs_menos_VSac_Tsubj,4); 
% freqStruct.VSpcs_menos_VSac.std = std(VSpcs_menos_VSac_Tsubj(~isnan(VSpcs_menos_VSac_Tsubj)),4); 

% linear regression
LinRegrTime.VSpcsNan.B1.mean = nanmean(B1,3);
LinRegrTime.VSpcsNan.B1.std = nanstd(B1,3);
LinRegrTime.VSpcsNan.B2.mean = nanmean(B2,3);
LinRegrTime.VSpcsNan.B2.std = nanstd(B2,3);
LinRegrTime.VSpcsNan.R2.mean = nanmean(R2,3);
LinRegrTime.VSpcsNan.R2.std = nanstd(R2,3);
LinRegrTime.VSpcsNan.p.mean = nanmean(p,3);
LinRegrTime.VSpcsNan.p.std = nanstd(p,3);

keyboard
fileOut   = fullfile(folderIn,'freqStruct500.mat');

save(fileOut,'freqStruct', '-v7.3')


%% freq vs fixRank
% mean amplitude in the time window of the fix duration 


%% mean across subjects
names = {'E01', 'E03', 'E04', 'E05', 'E07','E08', 'E09', 'E10',  'E11', 'E12', 'E14', 'J03', 'J04', 'J08', 'J11', 'J12'};

folderIn = [folderOut 'multipleFreq500/'];
freqStruct  = fullfile(folderIn,'freqStruct500.mat');
load(freqStruct);
keyboard

meanFreqStruct = structfun(@mean, freqStruct);

keyboard
chan = {'PO7','Fz','PO8','CPz','PO7','Pz','PO8','O1','Oz','O2','F1','F2','Cz'};

inFile = [folderIn 'EEGEYEE01_analysis-500_freq_3.set'];
EEG = pop_loadset(inFile);
times = EEG.times;

channels = {EEG.chanlocs.labels};
indchans = find(ismember({EEG.chanlocs.labels},chan));

keyboard

%for ch = 1:length(indchans)
for ch = 1:length(channels)   
    VSpcsAmpl = []; VSpcsNanAmpl = []; VSacAmpl = []; EXAmpl = []; VSac_EXAmpl = [];
    for min_freq = 3:30
        VSpcslala2 = []; VSpcsNanlala2 = []; VSaclala2 = []; EXlala2 = []; VSac_EXlala2 = [];
        for iSuj = 1:length(names);
            VSpcslala = freqStruct(iSuj).(strcat('VSpcs_freq',num2str(min_freq)));
            VSpcslala2 = [VSpcslala2; VSpcslala(ch,:)];
            
            VSpcsNanlala = freqStruct(iSuj).(strcat('VSpcsNan_freq',num2str(min_freq)));
            VSpcsNanlala2 = [VSpcsNanlala2; VSpcsNanlala(ch,:)];
            
            VSaclala = freqStruct(iSuj).(strcat('VSac_freq',num2str(min_freq)));
            VSaclala2 = [VSaclala2; VSaclala(ch,:)];

            EXlala = freqStruct(iSuj).(strcat('EX_freq',num2str(min_freq)));
            EXlala2 = [EXlala2; EXlala(ch,:)];
            
            VSac_EXlala = freqStruct(iSuj).(strcat('VSacmenosEX_freq',num2str(min_freq)));
            VSac_EXlala2 = [VSac_EXlala2; VSac_EXlala(ch,:)];
        end
        % mean across subjects
        VSacAmpl = [VSacAmpl; mean(VSaclala2,1)];
        EXAmpl = [EXAmpl; mean(EXlala2,1)];
        VSac_EXAmpl = [VSac_EXAmpl; mean(VSac_EXlala2,1)];
        
    end
    meanFreqStruct.(strcat('VSpcs_ch',char(channels(ch)))) = VSpcsAmpl;
    meanFreqStruct.(strcat('VSpcsNan_ch',char(channels(ch)))) = VSpcsNanAmpl;
    meanFreqStruct.(strcat('VSac_ch',char(channels(ch)))) = VSacAmpl;
    meanFreqStruct.(strcat('EX_ch',char(channels(ch)))) = EXAmpl;
    meanFreqStruct.(strcat('VSacmenosEX_ch',char(channels(ch)))) = VSac_EXAmpl;
end

fileOut   = fullfile(folderIn,'meanFreqStruct500VSacEX.mat');

save(fileOut,'meanFreqStruct')

%% time linear regression

names = {'E01', 'E03', 'E04', 'E05', 'E07','E08', 'E09', 'E10',  'E11', 'E12', 'E14', 'J03', 'J04', 'J08', 'J11', 'J12'};

folderIn = [folderOut 'multipleFreq500/'];
freqStruct  = fullfile(folderIn,'meanFreqStruct500VSacEX.mat');
load(freqStruct);

chan = {'PO7','Fz','PO8','CPz','PO7','Pz','PO8','O1','Oz','O2','F1','F2','Cz'};

inFile = [folderIn 'EEGEYEE01_analysis-500_freq_3.set'];
EEG = pop_loadset(inFile);
times = EEG.times;

indchans = find(ismember({EEG.chanlocs.labels},chan));
channels = {EEG.chanlocs.labels};


%% histo 2D

names = {'E01', 'E03', 'E04', 'E05', 'E07','E08', 'E09', 'E10',  'E11', 'E12', 'E14', 'J03', 'J04', 'J08', 'J11', 'J12'};

folderIn = [folderOut 'multipleFreq500/'];
freqStruct  = fullfile(folderIn,'meanFreqStruct500VSacEX.mat');
load(freqStruct);

chan = {'PO7','Fz','PO8','CPz','PO7','Pz','PO8','O1','Oz','O2','F1','F2','Cz'};

inFile = [folderIn 'EEGEYEE01_analysis-500_freq_3.set'];
EEG = pop_loadset(inFile);
times = EEG.times;

indchans = find(ismember({EEG.chanlocs.labels},chan));
channels = {EEG.chanlocs.labels};

for ch = 1:length(channels) 
    minVSac = min(meanFreqStruct.(strcat('VSac_ch',char(channels(ch))))(:));
    minEX = min(meanFreqStruct.(strcat('EX_ch',char(channels(ch))))(:));
    Min = min(minVSac,minEX);
    
    maxVSac = max(meanFreqStruct.(strcat('VSac_ch',char(channels(ch))))(:));
    maxEX = max(meanFreqStruct.(strcat('EX_ch',char(channels(ch))))(:));
    Max = max(maxVSac,maxEX);    
    
    figure(ch);
    imagesc(times,1:27,flipud(meanFreqStruct.(strcat('VSac_ch',char(channels(ch))))),[Min,Max]);colormap(bluewhitered(256)); colorbar
    xlabel('time (ms)')
    ylabel('freq (Hz)')
    title([channels(ch) '(VSac 500)'])
    colorbar
    % view(2)
    filename = strcat('VSac_chan',channels(ch),'.png');
    fileOut   = fullfile(folderIn,filename);
    y = 0:5:27;
    yTickLabels = arrayfun(@num2str,sort(y,'descend'),'uni',false);
    ax = gca;
    ax.YAxis.TickLabels = yTickLabels;
    saveas(gcf,fileOut{1})
end

for ch = 1:length(channels) 
    figure(ch);
    minVSac = min(meanFreqStruct.(strcat('VSac_ch',char(channels(ch))))(:));
    minEX = min(meanFreqStruct.(strcat('EX_ch',char(channels(ch))))(:));
    Min = min(minVSac,minEX);
    
    maxVSac = max(meanFreqStruct.(strcat('VSac_ch',char(channels(ch))))(:));
    maxEX = max(meanFreqStruct.(strcat('EX_ch',char(channels(ch))))(:));
    Max = max(maxVSac,maxEX);    
    
    imagesc(times,1:27,flipud(meanFreqStruct.(strcat('EX_ch',char(channels(ch))))),[Min,Max]);colormap(bluewhitered(256)); colorbar
    xlabel('time (ms)')
    ylabel('freq (Hz)')
    title([channels(ch) '(EX 500)'])
    colorbar
    % view(2)
    filename = strcat('EX_chan',channels(ch),'.png');
    fileOut   = fullfile(folderIn,filename);
    y = 0:5:27;
    yTickLabels = arrayfun(@num2str,sort(y,'descend'),'uni',false);
    ax = gca;
    ax.YAxis.TickLabels = yTickLabels;
    saveas(gcf,fileOut{1})
end

for ch = 1:length(channels) 
    figure(ch);
    minVSacEX = min(meanFreqStruct.(strcat('VSacmenosEX_ch',char(channels(ch))))(:));
    maxVSacEX = max(meanFreqStruct.(strcat('VSacmenosEX_ch',char(channels(ch))))(:));
    
    imagesc(times,1:27,flipud(meanFreqStruct.(strcat('VSacmenosEX_ch',char(channels(ch))))),[minVSacEX,maxVSacEX]);colormap(bluewhitered(256)); colorbar
    xlabel('time (ms)')
    ylabel('freq (Hz)')
    title([channels(ch) '(VS-EX 500)'])
    colorbar
    % view(2)
    filename = strcat('VSacmenosEX_chan',channels(ch),'.png');
    fileOut   = fullfile(folderIn,filename);
    y = 0:5:27;
    yTickLabels = arrayfun(@num2str,sort(y,'descend'),'uni',false);
    ax = gca;
    ax.YAxis.TickLabels = yTickLabels;
    saveas(gcf,fileOut{1})
end

%% Functions 

function EEGdata = fun_baselineCorrection(EEGdata,times)
    data_base = EEGdata(:,min(times)<times & times<-50,:);
    mean_data_base = mean(data_base,2);
    EEGdata = EEGdata - mean_data_base;
end

function especific_data = fun_freqChanTime(EEGdata,tabla,especific_data)
    lala_data = EEGdata(1:64,:,tabla.trial_number);
    lala = mean(lala_data,3);
    especific_data =  cat(3, especific_data, lala);
end

function especific_data = fun_freqChanTime_nan(EEGdata,tabla,especific_data,times)
    time2target = tabla.tonset;
    lala_data = EEGdata(1:64,:,tabla.trial_number);
    lala_data_nan = NaN(64,length(times),length(tabla.trial_number));
    for tr = 1:length(tabla.trial_number)
        ind = find([times]< 1000*time2target(tr));
        ind = ind(end);
        lala_data_nan(:,1:ind,tr) = lala_data(:,times<1000*time2target(tr));
    end
    lala = mean(lala_data_nan,3);
    especific_data =  cat(3, especific_data, lala);
end