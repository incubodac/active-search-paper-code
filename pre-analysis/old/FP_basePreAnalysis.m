classdef FP_basePreAnalysis
properties 
    Mode  % 'default', 'dac', 'bb'
    cfg   % configuration struct - Example under construction
    fn    % custom functions - each mode load from another FP_'mode'
end
methods
    % experimental constructor
    function obj   = FP_basePreAnalysis(varargin)
        % Constructor for the class
        % Mode: 'dac', 'bb', 'default'. different preprocessings
        % cfg: configuration struct. See default using only Mode argument
        if length(varargin) == 2
            obj.Mode = varargin{1};
            obj.cfg  = varargin{2};
        elseif length(varargin) == 1
            obj.Mode = varargin{1};
            disp('Loading default configurations')
            obj.cfg  = 'HACER';           
        end
        
        % Load custom functions for each mode
        if strcmpi(obj.Mode,'dac')
            fprintf('----------------------------------------------\n')
            fprintf('\t\trepo: FP_basePreAnalysis: prehybridana\n')
            fprintf('----------------------------------------------\n')
            obj.fn = FP_dac();
            obj.cfg = obj.fn.cfg;

        elseif strcmpi(obj.Mode,'bb')
            obj.fn = FP_bb();
            if length(varargin) == 1
                obj.cfg = obj.fn.cfg;
            end
        elseif strcmpi(obj.Mode,'default')
            print('Default mode, no custom functions')
        end

    end

    % preanalysis methods non static
    function ET         = ET_preprocessing(obj)
        % ET = ET_preprocessing()    
        % Eye tracking preProcessing. 
        % obj.cfg.fileIn  = ASCII file from eye tracker
        % obj.cfg.fileOut = complete file path for the mat output (output must be name as *.mat)
        % obj.cfg.keyword = keyword for parseeyelink to find events

        % mode :'default' and 'bb' use unmodified ET.events 
        % mode :'dac' modified hardcoded version marks are 290 for eyemap and 300 for end of

        fileIn  = obj.cfg.inFile;
        fileOut = obj.cfg.outFile; 
        keyword = obj.cfg.keyword;
        
        printLoadingFile(fileIn)
        
        if strcmp(obj.Mode,'default') || strcmp(obj.Mode,'bb')
            ET = parseeyelink(fileIn, fileOut, keyword);  
        elseif strcmp(obj.Mode,'dac') 
            ET = obj.fn.renameEtEvents(fileIn, fileOut, keyword);  
        end    
        
    end
    function EEG        = EEG_preprocessing(obj)
        % Load EDF file (pop_biosig) using ref vector as references 
        % Erase sensor channels and extra head channels 
        % Downsample to 1024Hz

        % obj.cfg.inFile = EDF file from EEG
        % obj.cfg.ref    = vector of references for pop_biosig
        % obj.cfg.nChansKeep = number of chans with data
        
        % mode = 'dac' rename channels other than first 64 channels
        % which are badly named in Nothingham experiment data
        % eg. :    EEG_Preprocessing(file,[65 66],'output_file.set', 64)
        
        inFile = obj.cfg.inFile;
        ref    = obj.cfg.ref;
        nChans = obj.cfg.nChansKeep;

        printLoadingFile(inFile)
        
        EEG = pop_biosig(inFile, 'ref', ref); % references are erased in the new file
        EEG = eeg_checkset(EEG);
        
        %channel names
        externalChan   = {'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8'};
        analogChan     = {'Ana1' 'Ana2' 'Ana3' 'Ana4' 'Ana5' 'Ana6' 'Ana7' 'Ana8'};
        sensorChan     = {'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp'};
        nbchan = EEG.nbchan; % number of channel in the data
        
        if strcmpi(obj.Mode,'dac')
            if nbchan > nChans      % renames external channels with wrong name
                for ch = nChans+1:nChans+8
                     EEG.chanlocs(ch).labels = externalChan{ch-nChans};
                end
             end
            if nbchan > nChans + 8  % renames analog channels with wrong name
                 for ch = nChans+9:nChans+16
                     EEG.chanlocs(ch).labels = analogChan{ch-(nChans+8)};
                 end
             end
            if nbchan > nChans      % erase extra channels
                EEG = pop_select(EEG,'nochannel',nChans+1:nbchan);
                EEG = eeg_checkset(EEG);
                fprintf('Unwanted channels erased\n')
            end
        end
        
        if strcmpi(obj.Mode,'dac') || strcmpi(obj.Mode,'bb')
            if EEG.srate > 1024 % Downsampling to 1024Hz
                EEG = pop_resample(EEG, 1024);
                EEG = eeg_checkset(EEG );
                fprintf('Resample to 1024hz\n')
            end
            if any(ismember({EEG.chanlocs.labels}, sensorChan))  % Delete sensor channels
                EEG = pop_select(EEG,'nochannel',sensorChan);
                EEG = eeg_checkset(EEG);
                fprintf('Erase sensors\n')
            end 
            if nbchan > 262 % Remove extra head channels
                % if nbchan > 262 I recorded more than 128 + external + analog
                chans256 = {EEG.chanlocs.labels}; % all channels
                chans128 = [chans256(1:128) chans256(end-15:end)]; % 128 + external + analog 
                EEG       = pop_select(EEG,'channel', chans128);
                EEG       = eeg_checkset(EEG);
                fprintf('Erase electrodes if save 256\n')
            end
        elseif strcmpi(obj.Mode,'default')
            fprintf('EDF with no modifications')
        end
        
        if nbchan < nChans  % Check if there are missing channels
            fprintf('Warning, subject with few channels\n')    
        end
        
        EEG = eeg_checkset(EEG);
    end
    function EEG        = EEG_prefiltering(obj)
        % Applies pass-band and notch filter  
        % notch filter is harcoded 47.5, 52.5, 846, 1,[],0
        % eg.: EEG_Preprocessing(file,[0.1 100],'output_file.set' , '/functions/resources/Standard-10-5-Cap385_witheog.elp')
        % Finally, loads chanlocs file to the EEG struct
        
        % obj.cfg.inFile = comple path to .set file
        % obj.cfg.preFilterEdges = edges for band pass filter using
        %                          pop_eegfiltnew()
        % obj.cfg.chanlocsFilePath = complete path to chanloc file
        % obj.cfg.noFiltAnalog = [Bool] wheter to remove analog before
        %                         filtering (1) or not (0)
        
        inFile           = obj.cfg.inFile;
        preFilterEdges   = obj.cfg.preFilterEdges;
        chanlocsFilePath = obj.cfg.chanlocsFilePath;
        noFiltAnalog     = obj.cfg.noFiltAnalog;
        
        % filtering power line noise. (output_file must be name as *.set)
        printLoadingFile(inFile)
        EEG                        = pop_loadset(inFile);
        EEG                        = eeg_checkset(EEG);  
        EEG.info.filter.locutoff   = preFilterEdges(1);
        EEG.info.filter.hicutoff   = preFilterEdges(2);
        
        if noFiltAnalog; analogData = EEG.data(end-7:end,:); end % remove analog data
        
        %notch filter
        EEG                        = pop_eegfiltnew(EEG, 47.5, 52.5, 846, 1,[],0);               
        
        %band pass filter
        EEG.info.filter.order      = 3*fix(EEG.srate/EEG.info.filter.locutoff); %DAC chequear origen de esto propuestas etc...
        EEG                        = pop_eegfiltnew(EEG, EEG.info.filter.locutoff, ...
                                                      EEG.info.filter.hicutoff, ...
                                                      EEG.info.filter.order, 0, [], 0);
        
        if noFiltAnalog; EEG.data = [EEG.data(1:end-8,:); analogData]; end %attaches unfiltered analog data
        
        EEG = eeg_checkset(EEG);
        
        if strcmpi(obj.Mode,'dac')
            % Add channel locations
            EEG = pop_chanedit(EEG, 'lookup', chanlocsFilePath);                                                  
            EEG = eeg_checkset(EEG);
        end
    end
    function EEG        = EEG_selectChanInterpolate(obj)
        % Plot EEG data by channel 
        % Plot channel with spectrum and variance of each channel
        % Use these plots to select which channels interpolate
        
        % obj.cfg.inFile: full path to .set EEG file
        % obj.cfg.nHeadChans: Number of channels to inspect before
        % interpolation
        
        inFile = obj.cfg.inFile;
        nChans = obj.cfg.nHeadChans;
        
        printLoadingFile(inFile)
        EEG = pop_loadset(inFile);
        EEG = eeg_checkset(EEG);

        % Delete analog for visual inspection
        chans = {EEG.chanlocs.labels};
        EEG   = pop_select(EEG,'channel', chans([1:nChans]));
        EEG   = eeg_checkset(EEG);

        eeglab redraw
        pop_eegplot( EEG, 1, 1, 1);
        
        figure; 
        subplot(1,3,1)
            pop_spectopo(EEG, 1, [], 'EEG' , 'freqrange',[0.1 100],'electrodes','off');
        subplot(1,3,2)
            %check by variance  
            chan=1:nChans;
            varch=[];
            for ch=1:nChans
                varch(ch) = var(EEG.data(ch,:));
            end
            plot( chan, varch,'*')
            set(gca, 'YScale','log')
            set(gca,'XTickLabel', [{EEG.chanlocs(1:20:128).labels}],'XTick',1:20:128)
            xlabel('Electrode')
            ylabel('log(variance)')
        subplot(1,3,3)
            plot(mean([EEG.data],2),'.')
            xticks(1:20:128)
            xticklabels([{EEG.chanlocs(1:20:128).labels}])
            xlabel('Electrode')
            ylabel('Mean Amplitud')

            
            
    end
    function EEG        = EEG_interpolate(obj)
        % Interpolate channels based on a list generated using information 
        % from EEG_selectChanInterpolate()

        % obj.cfg.inFile: full path to .set EEG file
        % obj.cfg.nHeadChans: Number of channels to interpolate
        % obj.cfg.badchanslist: cell with {'subj', [badchans]}
        % obj.cfg.CHANS: [MODE bb] Struct with chan information
        
        inFile = obj.cfg.inFile;
        badchanslist = obj.cfg.badchanslist;
        nChans = obj.cfg.nHeadChans;

        suj = strsplit(inFile,'/');
        suj = suj(end);
        suj = suj{1}(1:end-4);

        printLoadingFile(inFile)
        EEG = pop_loadset(inFile);
        EEG = eeg_checkset(EEG);
        badchans = badchanslist{strcmp({badchanslist{:,1}},suj),2};
        if ~isempty(badchans)
            fprintf('%5s ---> Interpolando...\n',suj)
            temp = pop_select( EEG,'nochannel',nChans+1:EEG.nbchan);
            temp = eeg_checkset( temp );
            
            if obj.cfg.nHeadChans == 128
                EEG = obj.fn.addChanInfo(EEG);
            end
            
            temp = eeg_interp(temp, badchans, 'invdist');
            temp = eeg_checkset( temp );

            EEG.data(badchans,:) = squeeze(temp.data(badchans,:));
        end
        EEG.info.badchans = badchans;
    end
    function EEG        = EEG_ET_sync(obj)
        % Load EEG and ET files and syn using pop_importeyetracker
        % obj.cfg.inFileEEG : full path to .set EEG file
        % obj.cfg.inFileET  : full path to .mat ET file
        % obj.cfg.marks     : marks to sync ET with EEG

        inFileEEG  = obj.cfg.inFileEEG;
        inFileET   = obj.cfg.inFileET;
        marks      = obj.cfg.marks;
        
        % Load EEG
        printLoadingFile(inFileEEG)
        EEG = pop_loadset(inFileEEG);
        EEG = eeg_checkset(EEG);  

        % Load ET
        printLoadingFile(inFileET)
        ET  = load(inFileET);

        if strcmpi(obj.Mode,'dac')
            fprintf('EEG events modified according to %s mode\n', obj.Mode)
            EEG = obj.fn.renameEEGevents(EEG);  
        elseif strcmpi(obj.Mode,'bb')
            if contains(inFileET,'pre')
                fprintf('EEG events modified according to %s mode\n', obj.Mode)
                [EEG, ET] = obj.fn.renameEEGevents(EEG, ET);
            end
        end  
        
        fprintf('\n')
        %fprintf('DAC (2022-05-06): select events: remove all events except 200 255\n')
        %EEG = pop_selectevent(EEG,'type',[200 255],'deleteevents','on');
        %fprintf('DAC (2022-05-06): resample EEG to 500Hz (same as ET)\n')
        EEG = pop_resample( EEG, 500);
        % synchronize eye-track
        fprintf('ESTO ES NUEVO\n')
        if EEG.srate ==512
            fprintf('DAC (2022-05-06): Synchro Radius = 4 samples (sr=512)\n')
            EEG = pop_importeyetracker(EEG, inFileET ,...
                                          marks,...
                                          1:size(ET.data,2),...
                                          ET.colheader,...
                                          1,1,1,1,4);
        elseif EEG.srate ==1024
            fprintf('DAC (2022-05-06): Synchro Radius = 8 samples (sr=1024)\n')
            EEG = pop_importeyetracker(EEG, inFileET ,...
                                          marks,...
                                          1:size(ET.data,2),...
                                          ET.colheader,...
                                          1,1,1,1,8);
        else
            fprintf('DAC (2022-05-06): Synchro Radius = 8 samples (sr=1024)\n')
            EEG = pop_importeyetracker(EEG, inFileET ,...
                                          marks,...
                                          1:size(ET.data,2),...
                                          ET.colheader,...
                                          1,1,1,1,5);
        end
        EEG = eeg_checkset(EEG);
        
        keyboard
    end
    function EEG        = engbertDetection(obj)
        % Bad data and eye movement detection using Engbert algorithms
        % The result of this function is Dimigen "rawEEG"
        
        % obj.cfg.inFile : full path to .set EEG file
        % obj.cfg.eye: eye to be used 'R' or 'L'
        % obj.cfg.screenSize: monitor resolution [x y]

        inFile  = obj.cfg.inFile;
        
        % Load EEG
        EEG = pop_loadset(inFile);
        EEG = eeg_checkset(EEG);

        if strcmp(obj.Mode,'bb')
            EEG = obj.fn.generateBestEyeChan(EEG);
        end

        eye        = obj.cfg.eye;
        screenSize = obj.cfg.screenSize;      
        
        gaze_x = find(strcmp({EEG.chanlocs.labels},[eye '-GAZE-X']));
        gaze_y = find(strcmp({EEG.chanlocs.labels},[eye '-GAZE-Y']));

        eyechans = [gaze_x gaze_y];

        EEG = pop_rej_eyecontin(EEG, eyechans, [1 1], ...
                                    screenSize, ...
                                    round(0.1*EEG.srate),...
                                    2); % 1: Reject bad data. 2: Mark bad data 
                                
        EEG = pop_detecteyemovements(EEG, [], eyechans, ...
                                        6, 4, 0.018364, ...
                                        1, 0, round(0.05*EEG.srate), 4, ...
                                        1, 1, 1);
        
        EEG = eeg_checkset(EEG);

    end
    function [sph, wts] = trainICA(obj)
        %if strcmpi(obj.Mode,'dac')
        % Filtrar a 2Hz - El prefiltro de arriba es [0.1 100] y notch. Para
        % ICA hacemos le filtramos <2Hz
        % 
        % obj.cfg.inFile: complete apth to the EEG file to load
        % obj.cfg.icaFilterEdge: edge for hi-pass filter
        % obj.cfg.trialDur: 
        % obj.cfg.preDur: 
        % obj.cfg.stimMark:  
        % obj.cfg.preMark:  
        % obj.cfg.nHeadChans
         
        inFile           = obj.cfg.inFile;
        icaFilterEdge    = obj.cfg.icaFilterEdge;
        sacLabel         = obj.cfg.sacLabel;

        nHeadChans       = obj.cfg.nHeadChans;
        conf              = obj.cfg;
        
        EEG             = obj.fn.prepareDataForICA(inFile, icaFilterEdge,conf);        
        badcahns         = EEG.info.badchans;

        %Dimingen sugested values
        OW_FACTOR        = 1;
        REMOVE_EPOCHMEAN = true;
        EEG_CHANNELS     = 1:nHeadChans; % indices of all EEG channels (exclude any eye-tracking channels here)
                                                 % Dimingen recommends to also include EOG channels (if also recorded against common reference)
                                                 %DAC: Since I don't have EOG channels I will use rescale
                                                 %eye-tracking data added
                                                 
        
        % DEFINIR como hacer reject bad data, doing this in epoched data
        % might result in great loss of data
        %-----------------------------------------------
        
        % overweight de SP
        EEG = pop_overweightevents(EEG, sacLabel,[-0.02 0.01],OW_FACTOR,REMOVE_EPOCHMEAN);         
        if strcmp(obj.Mode, 'dac')
            % TRAIN ICA (PCA?)
            if isempty(badcahns)
                fprintf('\nRunning ICA on optimized training data...')
                EEG = pop_runica(EEG,'extended',1,'interupt','on','chanind',EEG_CHANNELS); % or use binary ICA for more speed
            else
                fprintf('\nRunning ICA + PCA for %d on optimized training data...', obj.cfg.nChans-numel(badcahns))
                EEG = pop_runica(EEG,'extended',1,'interupt','on','chanind',EEG_CHANNELS,'options',{'pca',obj.cfg.nChans-numel(badcahns)}); % or use binary ICA for more speed
            end
        elseif strcmp(obj.Mode, 'bb')
            EEG = pop_runica(EEG, ...
                            'icatype','binica',...
                            'dataset',1,...
                            'chanind',EEG_CHANNELS,...
                            'options',{'extended',1,'pca',64});
        end
        
        % Remember ICA weights & sphering matrix 
        wts = EEG.icaweights;
        sph = EEG.icasphere;
    end
    function EEG        = epochData(obj)
        % Filtrar para analisis - Viene con [0.1 100], lo pasamos a [0.2 40]
         
        inFile           = obj.cfg.inFile;
        postFilterEdges  = obj.cfg.postFilterEdges;
        %chanlocsFilePath = obj.cfg.chanlocsFilePath;
        noFiltAnalog     = obj.cfg.noFiltAnalog;
        stimMark         = obj.cfg.stimMark; 
        preMark          = obj.cfg.preMark; 
        trialDur         = obj.cfg.trialDur;
        
        printLoadingFile(inFile)
        EEG                        = pop_loadset(inFile);
        %keyboard
        EEG                        = eeg_checkset(EEG);            
        EEG.info.filter.locutoff   = postFilterEdges(1);
        EEG.info.filter.hicutoff   = postFilterEdges(2);
        EEG.info.filter.order      = 3*fix(EEG.srate/EEG.info.filter.locutoff); %DAC chequear origen de esto propuestas etc...
        EEG                        = pop_eegfiltnew(EEG, EEG.info.filter.locutoff, ...
                                                      EEG.info.filter.hicutoff, ...
                                                      EEG.info.filter.order, 0, [], 0);
        EEG                        = eeg_checkset(EEG);                                                          
        % Epoch de estimulos - Inicio del trial [-0.2 TrialDur]
        %%%%%%%%%%%%%%DAC 2022 may 24 %%%%%%%%%%
        %%% We dont epoch now bc variable trial duration, we take eyemaps
        %%% out.
        %EEG                        = pop_epoch(EEG,{stimMark},[-0.2 trialDur]);
        
        event = EEG.event;
        event(~ismember({event.type},{'151','255'})) = [];
        Eind = find(strcmp({event.type},'151'));
        Sind = find(strcmp({event.type},'255'));
        eyemaps_start = [event(Eind([1:5:35])).latency]; 
        stims_start = [event(Sind([1:30:210])).latency];

        eyemaps_ranges = [eyemaps_start(1) stims_start(1)-100;...
                          eyemaps_start(2) stims_start(2)-100;...
                          eyemaps_start(3) stims_start(3)-100;...
                          eyemaps_start(4) stims_start(4)-100;...
                          eyemaps_start(5) stims_start(5)-100;...
                          eyemaps_start(6) stims_start(6)-100;...
                          eyemaps_start(7) stims_start(7)-100];
        %%%%%%%%%%%%%%%%%%%%%%
        %keyboard
        EEG = pop_select(EEG, 'nopoint', eyemaps_ranges,'channel', 1:148)
        
        
        EEG                        = eeg_checkset(EEG);                                                                 
    end
    function [EEGpre, EEGpos]   = epochDataPrePos(obj)
        % Filtrar para analisis - Viene con [0.1 100], lo pasamos a [0.2 40]
        inFile           = obj.cfg.inFile;
        postFilterEdges  = obj.cfg.postFilterEdges;
        chanlocsFilePath = obj.cfg.chanlocsFilePath;
        noFiltAnalog     = obj.cfg.noFiltAnalog;
        preTrialMark     = obj.cfg.preTrialMark; 
        posTrialMark     = obj.cfg.posTrialMark; 
        trialDur         = obj.cfg.trialDur;
 
        printLoadingFile(inFile{1})
        EEG                        = pop_loadset(inFile);
        EEG                        = eeg_checkset(EEG);            
        EEG.info.filter.locutoff   = postFilterEdges(1);
        EEG.info.filter.hicutoff   = postFilterEdges(2);
        EEG.info.filter.order      = 3*fix(EEG.srate/EEG.info.filter.locutoff); %DAC chequear origen de esto propuestas etc...
        EEG                        = pop_eegfiltnew(EEG, EEG.info.filter.locutoff, ...
                                                      EEG.info.filter.hicutoff, ...
                                                      EEG.info.filter.order, 0, [], 0);
        EEG                        = eeg_checkset(EEG);                                                          
        % Epoch de estimulos - Inicio del trial [-0.2 TrialDur]
        
        inele = ismember({EEG.event.type},'11');   
        intwe = ismember({EEG.event.type},'12');
        inthi = ismember({EEG.event.type},'13');
        infou = ismember({EEG.event.type},'14');
        infif = ismember({EEG.event.type},'15');
        inele = find(inele); 
        intwe = find(intwe); 
        inthi = find(inthi);
        infou = find(infou);
        infif = find(infif);
        goodind = [];
        badtwe = [];
        goodtwe = [];
        % |11|----1.5s---|12|--1.25s-1.75s----|13|----------4.5s-----------|14|--1.5s--|15|--1.5s---|16|           
            
        for trial_eleven = 1%1:numel(timeele)

            mask = timeele(trial_eleven)<[EEG.event.latency]/512 & [EEG.event.latency]/512 < (timeele(trial_eleven) + 10.75); 
            mask = find(mask);
            events = ismember({EEG.event(mask).type},{'12','13','14','15','16'});

            trialmarks = numel(unique({EEG.event(mask(events)).type}));

            if trialmarks==4  %trialmarks should be 4 to be considered a complete trail
                goodind 

            end

            end

        %             for twe_ev = intwe
%                 for fif_ev = infif
%                     time_distance2 = EEG.event(fif_ev).latency - EEG.event(ev).latency ;
%                     time_distance1 = EEG.event(ev).latency - EEG.event(twe_ev).latency ;
%                     if  time_distance1/512 < 1.76 && time_distance1/512 > 1.2 && time_distance2/512 < 5.55 && time_distance2/512 > 5.45
% 3
%                         goodind   = [goodind ev];
%                         goodtwe   = [goodtwe twe_ev];
%                         goodfif   = [goodfif fif_ev];
% 
% 
% 
%                     end
%                 end
%             end
        badtwe = setdiff(intwe,goodtwe);
        badfif = setdiff(infif,goodfif);
        EEG.event(badtwe) = [];
        EEG.event(badtwe) = [];
        badind = setdiff(inthi,goodind);
        badtrial = find(ismember(inthi,badind));
        EEG.info.badtrial = badtrial;

        EEGpre                        = pop_epoch(EEG,{preTrialMark},[-0.5 1.2]);
        EEGpre                        = eeg_checkset(EEGpre);
        EEGpos                        = pop_epoch(EEG,{posTrialMark},[-0.5 1.2]);
        EEGpos                        = eeg_checkset(EEGpos);    
        
        
    end
    function EEG        = filteringAnaData(obj)
        % Filtrar para analisis - Viene con [0.1 100], lo pasamos a [0.2 40]
         
        inFile           = obj.cfg.inFile;
        postFilterEdges  = obj.cfg.postFilterEdges;
        chanlocsFilePath = obj.cfg.chanlocsFilePath;
        noFiltAnalog     = obj.cfg.noFiltAnalog;
        stimMark         = obj.cfg.stimMark; 
        preMark          = obj.cfg.preMark; 
        trialDur         = obj.cfg.trialDur;
        
        %printLoadingFile(inFile{1})
        EEG                        = pop_loadset(inFile);
        EEG                        = eeg_checkset(EEG);            
        EEG.info.filter.locutoff   = postFilterEdges(1);
        EEG.info.filter.hicutoff   = postFilterEdges(2);
        EEG.info.filter.order      = 3*fix(EEG.srate/EEG.info.filter.locutoff); %DAC chequear origen de esto propuestas etc...
        EEG                        = pop_eegfiltnew(EEG, EEG.info.filter.locutoff, ...
                                                      EEG.info.filter.hicutoff, ...
                                                      EEG.info.filter.order, 0, [], 0);
        EEG                        = eeg_checkset(EEG);                                                          
        % Epoch de estimulos - Inicio del trial [-0.2 TrialDur]
        %EEG                        = pop_epoch(EEG,{stimMark},[-0.2 trialDur]);
        %EEG                        = eeg_checkset(EEG);                                                                 
    end

    function EEG        = epochDataPos(obj)
        % Filtrar para analisis - Viene con [0.1 100], lo pasamos a [0.2 40]
         
        inFile           = obj.cfg.inFile;
        postFilterEdges  = obj.cfg.postFilterEdges;
        chanlocsFilePath = obj.cfg.chanlocsFilePath;
        noFiltAnalog     = obj.cfg.noFiltAnalog;
        stimMark         = obj.cfg.stimMark; 
        posTrialMark     = obj.cfg.posTrialMark; 
        trialDur         = obj.cfg.trialDur;
        
        printLoadingFile(inFile{1})
        EEG                        = pop_loadset(inFile);
        EEG                        = eeg_checkset(EEG);            
        EEG.info.filter.locutoff   = postFilterEdges(1);
        EEG.info.filter.hicutoff   = postFilterEdges(2);
        EEG.info.filter.order      = 3*fix(EEG.srate/EEG.info.filter.locutoff); %DAC chequear origen de esto propuestas etc...
        EEG                        = pop_eegfiltnew(EEG, EEG.info.filter.locutoff, ...
                                                      EEG.info.filter.hicutoff, ...
                                                      EEG.info.filter.order, 0, [], 0);
        EEG                        = eeg_checkset(EEG);                                                          
        % Epoch de estimulos - Inicio del trial [-0.2 TrialDur]
        EEG                        = pop_epoch(EEG,{posTrialMark},[-0.5 1.2]);
        EEG                        = eeg_checkset(EEG);                                                                 
    end
    function EEG        = selectIcPlochl(obj)

        % Dimingen recomended parameters (hardcoded)
        IC_THRESHOLD = 1.1;   % variance ratio threshold (determined as suitable in Dimigen, 2019)
        SACC_WINDOW  = [2 0]; % saccade window (in samples!) to compute variance ratios (see Dimigen, 2019)
        PLOTFIG      = true;  % plot a figure visualizing influence of threshold setting?
%         ICPLOTMODE   = 1;     % plot component topographies (inverse weights)? (2 = only plot "bad" ocular ICs)
        FLAGMODE     = 3;     % overwrite existing rejection flags? (3 = yes)

        % Agarra el EEG de epochData y los datos de trainICA [spa wts]
        inFile    = obj.cfg.inFile;  % here it should be EEG file with cont. data and Engbert detection already done
        weights   = obj.cfg.weights; % struct saved after running trainICA() function          
        eegChans  = 1:obj.cfg.nHeadChans;
        preMark   = obj.cfg.preMark; 
        preDur    = obj.cfg.preDur;
        sacstring = obj.cfg.sacstring;
        fixstring = obj.cfg.fixstring;
        ICPLOTMODE = obj.cfg.plot;
        
        EEG = pop_loadset(inFile);
        EEG = eeg_checkset(EEG);      
        
        if strcmp(obj.Mode,'bb')
            EEG = obj.fn.preparePre(EEG, inFile);
            if obj.cfg.nHeadChans == 128
                EEG = obj.fn.addChanInfo(EEG);
            end
        end
        if strcmp(obj.Mode,'dac')  
           EEG = obj.fn.addChanInfo(EEG);
        end
        % Epoch pre experiments - Inicio del trial [-0.2 TrialDur]
        %EEG  = pop_epoch(EEG,{preMark},[-0.2 preDur]);
        EEG  = eeg_checkset(EEG);    

        % Remove any existing ICA solutions from your original dataset
        EEG = cleanICA(EEG);
        %DAC 2022 may 24
        EEG = pop_select(EEG,'channel',[1:128]);
        % load trained ica data within pre EEG struct
        EEG = loadICA(EEG, weights);

        %keyboard
        EEG.icachansind = eegChans; %!!!!!QUE PASA SI TENGO PCA CHEQUEAR DONDE SE USA ESTA VARIABLE
        EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv
        % Automatically flag ocular ICs (Plöchl et al., 2012)

        [EEG, ~] = pop_eyetrackerica(EEG, sacstring, fixstring,...
                                          SACC_WINDOW, IC_THRESHOLD,...
                                          FLAGMODE, PLOTFIG,...
                                          ICPLOTMODE);
    end
    function EEG        = selectIcManual(obj)
        eegChans = 1:obj.cfg.nHeadChans;
        weights  = obj.cfg.weights;  % struct saved after running trainICA() function          
        inFile   = obj.cfg.inFile;   % Feed with epoched data file
        badcomps = obj.cfg.badcomp;
        
        %load EEG epoched data using epochData() function
        EEG = pop_loadset(inFile); 
        
        if strcmp(obj.Mode,'bb')
            if max(eegChans) == 128
                EEG = obj.fn.addChanInfo(EEG);
            end
        end
        
        % Remove any existing ICA solutions from your original dataset
        EEG = cleanICA(EEG);

        % load trained ica data within pre EEG struct
        EEG = loadICA(EEG, weights);
        
        EEG.reject.gcompreject = badcomps;
        pop_selectcomps(EEG,1:size(EEG.icaweights,1)); 
    end
    function EEG        = removeComp(obj)
        %remove flagged ocular ICs
        eegChans = 1:obj.cfg.nHeadChans;
        weights  = obj.cfg.weights;  % struct saved after running trainICA() function          
        inFile   = obj.cfg.inFile;   % Feed with epoched data file
        badcomps = obj.cfg.badcomp;
        
        %load EEG epoched data using epochData() function
        EEG = pop_loadset(inFile); 
        
        
        
        if strcmp(obj.Mode,'bb')
            if max(eegChans) == 128
                EEG = obj.fn.addChanInfo(EEG);
            end
        end
        if strcmp(obj.Mode,'dac')  
             %keyboard
           extraChan = EEG.data(129:end,:);         %%%DAC 2022 may 24
           EEG = obj.fn.addChanInfo(EEG);
           EEG = pop_select(EEG,'channel',[1:128]); %%%DAC 2022 may 24
        end
        % Remove any existing ICA solutions from your original dataset
        EEG = cleanICA(EEG);

        % load trained ica data within pre EEG struct
        EEG = loadICA(EEG, weights);
        
        %keyboard
        EEG.icachansind = eegChans; %!!!!!QUE PASA SI TENGO PCA CHEQUEAR DONDE SE USA ESTA VARIABLE
        EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv
        
        EEG = pop_subcomp(EEG,find(badcomps)); % remove ICs
        %%%DAC 2022 may 24
        EEG.data = cat(1,EEG.data,extraChan);
        EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv
    end



    function qualityCheck(obj)
        % Agarra el EEG de epochData y los datos de trainICA [spa wts]
        
        % EEG.icasphere = spa 
        % EEG.icaweigths = wts 
    end
    function EEG       = export_events2csv(obj)

        ETfile = obj.cfg.ETmat;
        out=regexp(ETfile,'/','split');
        fprintf('Loading ET messages from subject %s', out{end})
        ET = load(ETfile);
        fix=[];
        vs_start =[];
        vs_end =[];
        mem_start =[];
        mem_end =[];


        for i = 1:numel(ET.event(:,1))
            if ET.event(i,2)==50
                fix =[fix ; ET.event(i,1)];
            elseif ET.event(i,2)==200
                mem_start =[mem_start ; ET.event(i,1)];
            elseif ET.event(i,2)==201
                mem_end =[mem_end ; ET.event(i,1)];
            elseif ET.event(i,2)==250
                vs_start =[vs_start ; ET.event(i,1)];
            elseif ET.event(i,2)==251
                vs_end =[vs_end ; ET.event(i,1)];
            end
        end
        vs_phase_dur  =  vs_end    - vs_start;


        inFile   = obj.cfg.inFile;   % Feed with 4analyzed data 
        out=regexp(inFile,'/','split');
        fprintf('Loading EEG analysis data from subject %s', out{end})
        EEG = pop_loadset(inFile); 
        
        evt_ind = length(EEG.event); %last event idx to add starting from there
        types = {'cross1','mem','cross2','vs'}; %list of events to add, 'vs' It's already there but needs to be renamed
        evt_vec = [vs_start - fix vs_start - mem_start vs_start - mem_end];% event times taken from vs_start trigger
        evt_dur_vec = [mem_start - fix mem_end - mem_start vs_start - mem_end];% event times taken from vs_start trigger

        eeg_vs_evt_idx = find(ismember({EEG.event.type},'250')); %EEG idx for vs_start

        i=0;
        for itr = eeg_vs_evt_idx
            i = i+1;
            vs_lat = int64(EEG.event(itr).latency);
            vs_end_time = EEG.times(vs_lat) + vs_phase_dur(i);
            vs_end_idx = nearestpoint(vs_end_time,EEG.times);
            duration =vs_end_idx-vs_lat;
            %change trigger from '250' to 'vs'
            EEG.event(itr).type = 'vs';
            EEG.event(itr).duration = double(duration);
            %loop to add events for other phases of the experiment
            for intr = 1:3
                evt_ind=evt_ind+1;%starting from last event in EEG.event
                EEG.event(evt_ind).type = types{intr};   
                evt_time = EEG.times(vs_lat) - evt_vec(i,intr);
                evt_lat = nearestpoint(evt_time,EEG.times);
                latency = evt_lat;
                duration = int64(evt_dur_vec(i,intr)*EEG.srate/1000);%vs_lat - evt_lat; may have 1ms error didnt use nearestpoint
                EEG.event(evt_ind).latency = latency; 
                EEG.event(evt_ind).duration = double(duration);
            end
        end
        EEG=eeg_checkset(EEG,'eventconsistency');%reorder events according to latency

        ev_et_times = [fix ;mem_start ;mem_end ;vs_start];
        ev_et_times = sort(ev_et_times);
        ind = find(ismember({EEG.event.type}, types));
        times = EEG.times(int64([EEG.event(ind).latency]));
        [cor ,p] = corr(times',ev_et_times);
        plot(times',ev_et_times,'.')      
        % Add legend with correlation and p-value
        legend(sprintf('Correlation = %.2f, p-value = %.2f', cor, p));
        title('Comparing experimental phases ET marks and EEG imported ones')
        
    end
    function            EEG = reattachETchannels(obj,EEG)
        inFile   =    obj.cfg.inFile;
        EEG4et   = pop_loadset(inFile);
        etchans  = {'R-GAZE-X' 'R-GAZE-Y' 'R-AREA'};
        %keyboard
        EEG.chanlocs = EEG4et.chanlocs;
        ind4newchan = find(ismember({EEG4et.chanlocs.labels},etchans));
        EEG.data(ind4newchan,:) = EEG4et.data(ismember({EEG4et.chanlocs.labels},etchans),:);     

    end
    
end

%static methods don't require to intatiate the class
methods(Static)
end %methods


end % classdef

%% Functions Definition
function printLoadingFile(fullPath)
    % check version name = split(EDF_file,'/');
    name = strsplit(fullPath,'/');
    fprintf(['Loading file: ' name{end} '\n'] )
end
function EEG = cleanICA(EEG)
    EEG.icaact      = [];
    EEG.icasphere   = [];
    EEG.icaweights  = [];
    EEG.icachansind = [];
    EEG.icawinv     = [];
    EEG = eeg_checkset(EEG); 
end
function EEG = loadICA(EEG, w)
        W = load(w);
        EEG.icasphere   = W.sph; 
        EEG.icaweights  = W.wts;
        EEG = eeg_checkset(EEG); 
end