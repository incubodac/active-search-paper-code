% Edited by jek 10/11/22 to add a new user

% function [runpath,code_path,session_path,whoisrunning, cfg]=add_paths_EEG_ET_beta(print_info,laptop)
%add_paths_matlab(print_info)  
% if(nargin == 0)
%     print_info = 1;
%     laptop = [];
% end
% if(nargin == 1)
%     laptop = [];
% end
% [~,whoisrunning] = system('whoami')

function [runpath,code_path,session_path,whoisrunning, cfg]=add_paths_EEG_ET_beta(print_info,whoisrunning)

if(nargin == 0)
    print_info = 1;
    [~,whoisrunning] = system('whoami');
end
if(nargin == 1)
    [~,whoisrunning] = system('whoami');
end

switch strtrim(whoisrunning);

     case 'D' %'dac-Desktop'
        if 1
            dropbox_path            = '/media/dac/data/repos/visualsearch_eegem/codes';
            local_path              = '/Volumes/DAC1T/Analysis_2020';
            local_out               = '/Volumes/DAC1T/Analysis_2020'
            code_path.eeglabpath    = '~/Documents/MATLAB/eeglab2023.0';
            code_path.psyco         ='/home/dac/Toolbox/psyco';
            code_path.heatmap       ='/media/dac/data/repos/visualsearch_eegem/codes/my_functions/heatmap_code'
            code_path.my_functions  = '/media/dac/data/repos/visualsearch_eegem/codes/my_functions'
            runpath                 = fullfile(dropbox_path,'');
            % to use as datapath in Bruno script
            code_path.raw           = '/media/dac/data/Experimentos/Analysis_2020/raw_BDF_data/'
            code_path.corregistro_fn = '/media/dac/data/repos/corregistro/codes'
            code_path.toolbox             = '/usr/local/MATLAB/R2016a/toolbox/'
            code_path.functions     = '';
  
        else
            if (laptop==1)
                dropbox_path            = '/Users/dac/Documents/GitLab/corregistro/Analysis_dac' %'/Volumes/DAC-DRIVE/codes';
            if 7== exist('/Volumes/DAC1T/Analysis_2020')
                local_path              = '/Volumes/DAC1T/Analysis_2020'
                local_out               = '/Volumes/DAC1T/Analysis_2020'
            elseif 7== exist('/Volumes/DAC1Tsolid/Analysis_2020')
                local_path              = '/Volumes/DAC1Tsolid/Analysis_2020'
                local_out               = '/Volumes/DAC1Tsolid/Analysis_2020'
            elseif (7==exist('/media/dac/DATOS-curie/')) & (0 == exist('/Volumes/DAC1Tsolid/Analysis_2020'))
                local_path              = '/media/dac/DATOS-curie/Analysis_2020'
                local_out               = '/media/dac/DATOS-curie/Analysis_2020'
            else 
                local_path              = '/Volumes/homes/admin/Analysis_2020'%'/Volumes/DAC1T/Analysis_2020'%'/Volumes/DAC-DRIVE/EEGEYE/Analysis_2020' old disk% /Volumes/homes/admin/Analysis_2020
                local_out               = '/Volumes/homes/admin/Analysis_2020'%'/Volumes/DAC1T/Analysis_2020'% old disk%'/Volumes/homes/admin/Analysis_2020'
            end
                code_path.eeglabpath    = '/Users/dac/Documents/MATLAB/eeglab14_1_2b'
                code_path.psyco         = '/Users/dac/Documents/MATLAB/psyco'
                code_path.heatmap       = '/Users/dac/Documents/Gitlab/visualsearch_eegem/codes/my_functions/heatmap_code'
            % to use as datapath in Bruno script

                code_path.raw           = [  local_path '/raw_BDF_data/']
                code_path.corregistro_fn = '/Users/dac/Documents/Gitlab/corregistro/codes'
                code_path.toolbox             = '/Users/dac/Documents/MATLAB'
                code_path.my_functions  = '/Users/dac/Documents/Gitlab/visualsearch_eegem/codes/my_functions'      
                code_path.functions  = '/Users/dac/Documents/Gitlab/corregistro/Analysis_dac/functions' 
                

            end
                    
            runpath                 = fullfile(dropbox_path,'');
        end
        
        % session paths  
        session_path.raw            = fullfile(local_path,'raw_BDF_data');
        session_path.rawet          = fullfile(local_path,'raw_ASC_data');
        session_path.matfiles       = fullfile(local_path,'mat_files/');
        session_path.renamed        = fullfile(local_path,'1.renamed/');
        session_path.prefiltered    = fullfile(local_out,'2.prefiltered/');
        session_path.out            = local_out
        session_path.interpolated    = fullfile(local_out,'3.interpolated/');
        session_path.et_mat         = fullfile(local_out,'ET_mat_files/');
        session_path.synced         = fullfile(local_out,'4.EEG_ET_synced/');
        session_path.eng            = fullfile(local_out,'5.raw_eng/');
        session_path.eog            = fullfile(local_out,'5.1raw_eng/');
        session_path.icaWeights     = fullfile(local_out,'6.ica_wts/'); 
        session_path.epoched        = fullfile(local_out,'7.epoched_data/');
        session_path.data_analysis  = fullfile(local_out,'8.data_analysis/');
        session_path.aux_data       = fullfile(local_out,'eprs/');
        session_path.conditions_data       = fullfile(local_out,'A1.conditions_epochs/');
        session_path.lmm            = fullfile(local_out, 'LMM/');

        session_path.subjname   = {'E01','E03','E04','E05','E06','E07','E08','E09',...
            'E10','E11','E12','E13','E14','J01','J02','J03','J04','J05','J06','J07','J08','J09','J10','J11','J12'};
        session_path.sessionfilenames = {'E010620';'E030622';'E040625';'E050625';'E060626';...
            'E070627';'E080628';'E090702';'E100703';...
            'E110703';'E120703';'E130705';'E140709';'J010705';'J020705';'J030705';'J040712';...
            'J050712';'J060714';'J070714';'J080714';'J090714';'J100715';'J110715';'J120715'};
        session_path.Trialfilenames = {'S10620_Trials';'E010622_Trials';'E040625_Trials';'E050625_Trials';'E060626_Trials';...
            'E070627_Trials';'E080628_Trials';'E080702_Trials';'E100703_Trials';...
            'E110703_Trials';'E120703_Trials';'E13missing';'E140709';'J120715_Trials'};
        session_path.whichEyeAnalyze     = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};
        session_path.whichEye            = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};
    % mdf 22/02/21: agrego un case con mis paths   
    case 'root'
        if laptop==2
            dropbox_path            = '/home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac' %'/Volumes/DAC-DRIVE/codes';
            local_path              = '/media/cbclab/MARIADAFON1T/Analysis_2020'% idealmente iria al NAS o a mi disco interno
            local_out               = '/media/cbclab/MARIADAFON1T/Analysis_2020'% idealmente iria al NAS o a mi disco interno
            code_path.eeglabpath    = '/usr/local/MATLAB/R2018b/toolbox/eeglab2021.0'
            code_path.psyco         = '/usr/local/MATLAB/R2018b/toolbox/psyco'
            code_path.heatmap       = '/home/cbclab/Documents/gitlab_liaa/visualsearch_eegem/codes/my_functions/heatmap_code'
        
            % to use as datapath in Bruno script
            code_path.raw           = [  local_path '/raw_BDF_data/']
            code_path.corregistro_fn = '/home/cbclab/Documents/gitlab_liaa/corregistro/codes'
            code_path.toolbox       = '/usr/local/MATLAB/R2018b/toolbox'
            code_path.my_functions  = '/home/cbclab/Documents/gitlab_liaa/visualsearch_eegem/codes/my_functions'         
            code_path.functions     = '/home/cbclab/Documents/gitlab_liaa/corregistro/Analysis_dac/functions' % mdf 2/3/21 para que FP_epochAnalysis.m esté acá
 
            runpath                 = fullfile(dropbox_path,'');
        
            %session paths    
            session_path.raw            = fullfile(local_path,'raw_BDF_data');
            session_path.rawet          = fullfile(local_path,'raw_ASC_data');
            session_path.matfiles       = fullfile(local_path,'mat_files/');
            session_path.renamed        = fullfile(local_path,'1.renamed/');
            session_path.prefiltered    = fullfile(local_out,'2.prefiltered/');
            session_path.out            = local_out
            session_path.interpolated    = fullfile(local_out,'3.interpolated/');
            session_path.et_mat         = fullfile(local_out,'ET_mat_files/');
            session_path.synced         = fullfile(local_out,'4.EEG_ET_synced/'); % mdf para plot_eyeMovs.m 
            session_path.eng            = fullfile(local_out,'5.raw_eng/');
            session_path.eog            = fullfile(local_out,'5.1raw_eng/');
            session_path.icaWeights     = fullfile(local_out,'6.ica_wts/'); 
            session_path.epoched        = fullfile(local_out,'7.epoched_data/');
            session_path.data_analysis  = fullfile(local_out,'8.data_analysis/');
            session_path.data_analysis500  = fullfile(local_out,'8.data_500/');
            
            % session_path.remove_badIC  =
            % fullfile(local_out,'9.removed_badIC/');% mdf 26/4/21 since we
            % couldnt rescue J10 participant, this step is unnecessary
            session_path.aux_data       = fullfile(local_out,'eprs/');
            session_path.conditions_data       = fullfile(local_out,'M_epochs_conditions/');
            session_path.tfce           = fullfile(local_out,'TFCE_results/');
            session_path.fixs           = fullfile(local_out,'fixs/');
            session_path.plots          = fullfile(local_out,'plots/');
            session_path.freqEEG        = fullfile(local_out,'freqEEG/'); % for oscillation analysis
            session_path.BrunoData      = fullfile(local_out,'brunobian/continous_data/'); 
            
            session_path.subjname   = {'E01','E03','E04','E05','E06','E07','E08','E09',...
                'E10','E11','E12','E13','E14','J01','J02','J03','J04','J05','J06','J07','J08','J09','J10','J11','J12'};
            session_path.sessionfilenames = {'E010620';'E030622';'E040625';'E050625';'E060626';...
                'E070627';'E080628';'E090702';'E100703';...
                'E110703';'E120703';'E130705';'E140709';'J010705';'J020705';'J030705';'J040712';...
                'J050712';'J060714';'J070714';'J080714';'J090714';'J100715';'J110715';'J120715'};
            session_path.Trialfilenames = {'S10620_Trials';'E010622_Trials';'E040625_Trials';'E050625_Trials';'E060626_Trials';...
                'E070627_Trials';'E080628_Trials';'E080702_Trials';'E100703_Trials';...
                'E110703_Trials';'E120703_Trials';'E13missing';'E140709';'J120715_Trials'};
            session_path.whichEyeAnalyze     = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};
            session_path.whichEye            = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};
        elseif laptop==3
            dropbox_path            = '/home/maria/Documents/gitlab_liaa/corregistro/Analysis_dac' %'/Volumes/DAC-DRIVE/codes';
            local_path              = '/media/maria/MARIADAFON1T/Analysis_2020'% idealmente iria al NAS o a mi disco interno
            local_out               = '/media/maria/MARIADAFON1T/Analysis_2020'% idealmente iria al NAS o a mi disco interno
            code_path.eeglabpath    = '/usr/local/MATLAB/R2018b/toolbox/eeglab2021.0'
            code_path.psyco         = '/usr/local/MATLAB/R2018b/toolbox/psyco'
            code_path.heatmap       = '/home/maria/Documents/gitlab_liaa/visualsearch_eegem/codes/my_functions/heatmap_code'

            % to use as datapath in Bruno script
            code_path.raw           = [  local_path '/raw_BDF_data/']
            code_path.corregistro_fn = '/home/maria/Documents/gitlab_liaa/corregistro/codes'
            code_path.toolbox       = '/usr/local/MATLAB/R2018b/toolbox'
            code_path.my_functions  = '/home/maria/Documents/gitlab_liaa/visualsearch_eegem/codes/my_functions'         
            code_path.functions     = '/home/maria/Documents/gitlab_liaa/corregistro/Analysis_dac/functions' % mdf 2/3/21 para que FP_epochAnalysis.m esté acá
 
            runpath                 = fullfile(dropbox_path,'');
            
            %session paths    
            session_path.raw            = fullfile(local_path,'raw_BDF_data');
            session_path.rawet          = fullfile(local_path,'raw_ASC_data');
            session_path.matfiles       = fullfile(local_path,'mat_files/');
            session_path.renamed        = fullfile(local_path,'1.renamed/');
            session_path.prefiltered    = fullfile(local_out,'2.prefiltered/');
            session_path.out            = local_out
            session_path.interpolated   = fullfile(local_out,'3.interpolated/');
            session_path.et_mat         = fullfile(local_out,'ET_mat_files/');
            session_path.synced         = fullfile(local_out,'4.EEG_ET_synced/'); % mdf para plot_eyeMovs.m 
            session_path.eng            = fullfile(local_out,'5.raw_eng/');
            session_path.eog            = fullfile(local_out,'5.1raw_eng/');
            session_path.icaWeights     = fullfile(local_out,'6.ica_wts/'); 
            session_path.epoched        = fullfile(local_out,'7.epoched_data/');
            session_path.data_analysis  = fullfile(local_out,'8.data_analysis/');
            session_path.data_analysis500  = fullfile(local_out,'8.data_500/');
            
            % session_path.remove_badIC  =
            % fullfile(local_out,'9.removed_badIC/');% mdf 26/4/21 since we
            % couldnt rescue J10 participant, this step is unnecessary
            session_path.aux_data       = fullfile(local_out,'eprs/');
            session_path.conditions_data       = fullfile(local_out,'M_epochs_conditions/');
            session_path.tfce           = fullfile(local_out,'TFCE_results/');
            session_path.fixs           = fullfile(local_out,'fixs/');
            session_path.plots          = fullfile(local_out,'plots/');
            session_path.freqEEG        = fullfile(local_out,'freqEEG/'); % for oscillation analysis
            session_path.BrunoData      = fullfile(local_out,'brunobian/continous_data/'); 
            
            session_path.subjname   = {'E01','E03','E04','E05','E06','E07','E08','E09',...
                'E10','E11','E12','E13','E14','J01','J02','J03','J04','J05','J06','J07','J08','J09','J10','J11','J12'};
            session_path.sessionfilenames = {'E010620';'E030622';'E040625';'E050625';'E060626';...
                'E070627';'E080628';'E090702';'E100703';...
                'E110703';'E120703';'E130705';'E140709';'J010705';'J020705';'J030705';'J040712';...
                'J050712';'J060714';'J070714';'J080714';'J090714';'J100715';'J110715';'J120715'};
            session_path.Trialfilenames = {'S10620_Trials';'E010622_Trials';'E040625_Trials';'E050625_Trials';'E060626_Trials';...
                'E070627_Trials';'E080628_Trials';'E080702_Trials';'E100703_Trials';...
                'E110703_Trials';'E120703_Trials';'E13missing';'E140709';'J120715_Trials'};
            session_path.whichEyeAnalyze     = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};
            session_path.whichEye            = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};
        end

    % jek 10/11/22: agrego un case con mis paths   
    case 'J'
        dropbox_path            = '/home/juank/Desktop/EJN2022';            %'/Volumes/DAC-DRIVE/codes';
        local_path              = '/home/juank/Desktop/EJN2022/matlab';     % idealmente iria al NAS o a mi disco interno
        local_out               = '/home/juank/Desktop/EEGEYE_EJN2022';     % idealmente iria al NAS o a mi disco interno
        code_path.eeglabpath    = '/home/juank/toolbox/eeglab2021.0';
        code_path.psyco         = '';%'/usr/local/MATLAB/R2018b/toolbox/psyco'
        code_path.heatmap       = '/home/juank/repos/visualsearch_eegem/codes/my_functions/heatmap_code';
        
        code_path.raw           = [  local_path '/raw_BDF_data/'];
        code_path.corregistro_fn= '/home/juank/repos/corregistro/codes';
        code_path.toolbox       = '/usr/local/MATLAB/R2018b/toolbox';
        code_path.my_functions  = '';%'/home/juank/repos/visualsearch_eegem/codes/my_functions'         
        code_path.functions     = '';%'/home/juank/repos/corregistro/Analysis_dac/functions' % mdf 2/3/21 para que FP_epochAnalysis.m esté acá

        runpath                 = fullfile(dropbox_path,'');
        
        %session paths    
        session_path.raw            = fullfile(local_path,'raw_BDF_data');
        session_path.rawet          = fullfile(local_path,'raw_ASC_data');
        session_path.matfiles       = fullfile(local_path,'mat_files/');
        session_path.renamed        = fullfile(local_path,'1.renamed/');
        session_path.prefiltered    = fullfile(local_out,'2.prefiltered/');
        session_path.out            = local_out;
        session_path.interpolated    = fullfile(local_out,'3.interpolated/');
        session_path.et_mat         = fullfile(local_out,'ET_mat_files/');
        session_path.synced         = fullfile(local_out,'4.EEG_ET_synced/'); % mdf para plot_eyeMovs.m 
        session_path.eng            = fullfile(local_out,'5.raw_eng/');
        session_path.eog            = fullfile(local_out,'5.1raw_eng/');
        session_path.icaWeights     = fullfile(local_out,'6.ica_wts/'); 
        session_path.epoched        = fullfile(local_out,'7.epoched_data/');
        session_path.data_analysis  = fullfile(local_out,'8.data_analysis/');
        session_path.data_analysis500  = fullfile(local_out,'8.data_500/');
        
        % session_path.remove_badIC  =
        % fullfile(local_out,'9.removed_badIC/');% mdf 26/4/21 since we
        % couldnt rescue J10 participant, this step is unnecessary
        session_path.aux_data       = fullfile(local_out,'eprs/');
        session_path.conditions_data= fullfile(local_out,'M_epochs_conditions/');
        session_path.tfce           = fullfile(local_out,'TFCE_results/');
        session_path.fixs           = '/home/juank/Desktop/EJN2022/data/output/'; %fullfile(local_out,'fixs/'); jek 10/11/2022: This should be automatic in Python...
        session_path.plots          = fullfile(local_out,'plots/');
        session_path.freqEEG        = fullfile(local_out,'freqEEG/'); % for oscillation analysis
        session_path.BrunoData      = fullfile(local_out,'brunobian/continous_data/'); 

        session_path.subjname   = {'E01','E03','E04','E05','E06','E07','E08','E09',...
            'E10','E11','E12','E13','E14','J01','J02','J03','J04','J05','J06','J07','J08','J09','J10','J11','J12'};
        session_path.sessionfilenames = {'E010620';'E030622';'E040625';'E050625';'E060626';...
            'E070627';'E080628';'E090702';'E100703';...
            'E110703';'E120703';'E130705';'E140709';'J010705';'J020705';'J030705';'J040712';...
            'J050712';'J060714';'J070714';'J080714';'J090714';'J100715';'J110715';'J120715'};
        session_path.Trialfilenames = {'S10620_Trials';'E010622_Trials';'E040625_Trials';'E050625_Trials';'E060626_Trials';...
            'E070627_Trials';'E080628_Trials';'E080702_Trials';'E100703_Trials';...
            'E110703_Trials';'E120703_Trials';'E13missing';'E140709';'J120715_Trials'};
        session_path.whichEyeAnalyze     = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};
        session_path.whichEye            = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right'};

    otherwise
        warning(['User unknown: ' whoisrunning ' , check paths']);
end

%       CONFIGURATION STRUCT TO USE WITH FP_basePreAnalysis

%general use for passinginput and ouput files to pre-processing fucntions,
cfg.inFile               = [];
cgf.outfile              = [];

%parameters
cfg.eye                  = 'R';                                 % eye to be used 'R' or 'L'
cfg.ref                  = [69 70];                             % vector of references for pop_biosig
cfg.nChans               = 64;                                  % number of chans with data
cfg.screenSize           = [1920 1080];                         % monitor resolution [x y]
cfg.chanlocsFilePath     = [code_path.eeglabpath ...
    '/functions/resources/Standard-10-5-Cap385_witheog.elp'];   % complete path to chanloc file               
%ET
cfg.keyword              = 'BUTTON';                            % keyword to parse eyetracker data

%EEG prefiltering
cfg.lineFreq             = 50;
cfg.filterEdges          = [0.1 100];                           % edges for band pass filter using pop_eegfiltnew()

cfg.noFiltAnalog         = 0;                                   % [Bool] wheter to remove analog before filtering (1) or not (0)

%interpolation
cfg.nHeadChans           = 64;                                  % Number of channels to inspect  before interpolation
cfg.badchanslist         = [];                                  % cell with {'subj', [badchans]} update  after inspetion

% Load EEG and ET files and syn using pop_importeyetracker
cfg.inFileEEG           = [];                                   % full path to .set EEG file
cfg.inFileET            = [];                                   % full path to .mat ET file
cfg.marks               = [290 300];                            % marks to sync ET with EEG

if(print_info)
    disp(whoisrunning)
    disp(code_path)
    disp(runpath)
end

end

