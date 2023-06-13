classdef FP_bb
properties
    cfg
end
methods
    % experimental constructor
    function obj = FP_bb()
        % Constructor for the class
        obj.cfg = struct();
    end
    function EEG = assignFix(obj)
        EEG = pop_loadset(obj.cfg.eegFile);
        EEG = eeg_checkset(EEG);
        ET = load(obj.cfg.etFile);

        % Load DATA and ESTIM
        load(obj.cfg.matFile);
        load('/home/brunobian/Documents/Repos/expe_corregistro_2018/Analisis_ET/files/ESTIM')

        fields = {'length' 'catgram' 'SntType' 'pos' 'wrdID' 'palabras'...
                   'paltarget' 'freq' 'pred' 'predPrev' 'predNext' ...
                   'RP' 'posRelRP'};

        % Add the estimuli information to the EEG DATA       
        DATA   = addFieldsEstim(DATA, ESTIM, fields);
        % Add position of spaces information 
        DATA   = addSpaces(DATA);
        % Add latancies (timepoints of senteces starts and ends)
        DATA   = addLatency(DATA, EEG);

        EEG.info.DATA = DATA;
        EEG.info.init = init;

        EEG = eeg_checkset(EEG);

        yLim = init.height*0.15;
        
        events  = {EEG.event.type};
        ends    = find(strcmp(events, '221'));
        fixX = [EEG.event.fix_avgpos_x];
        fixY = [EEG.event.fix_avgpos_y] - init.height/2;
        
        
        clear ET DATA init DATAPRACT ESTIM
        for i = 1:length(EEG.info.DATA)

            % Find cross for this trial
            cross    = find(strcmp(events, num2str(i)),1,'last');
            % Find all text starts (230: proverbs - 231: common)
            redPoint = find(strcmp(events, '230') + strcmp(events, '231'));
            % Find text start for this trial (first after cross)
            trStart = redPoint(find((redPoint - cross) > 0,1));
            trEnd   = ends(find((ends - trStart)>0,1));

            fixTr = find(strcmp(events, ['B_fixation']));
            fixTr = fixTr(fixTr > trStart & fixTr < trEnd);

            xPos = fixX(fixTr);
            yPos = fixY(fixTr);
            edges = EEG.info.DATA(i).spaces;

            [~, bin] = histc(xPos, edges);

            % Delete fixation outside any word and transform in index (0
            % starting)
            
            badFix = bin>0 & abs(yPos)<yLim;
            fixTr = fixTr(badFix);
            bin = bin(badFix)-1; 

            EEG.info.DATA(i).fixAsigWrd = bin;         
            EEG.info.DATA(i).fixAsigID  = EEG.info.DATA(i).ind * 1000 + bin;
            EEG.info.DATA(i).fixEventInd = fixTr;
            EEG.info.DATA(i).fixRank     = 1:length(bin);

            if isempty(bin)
                suj = obj.cfg.eegFile(end-6:end-4);
                fprintf('NO FIX: Trial %.0f, Subject %s\n', i, suj)
                continue
            end

            % Generate First Pass Reading Fixations (FPRF)
            indGrow = logical([1 diff(bin)>=0]);

            EEG.info.DATA(i).FprfAsigW   = EEG.info.DATA(i).fixAsigWrd(indGrow); 
            EEG.info.DATA(i).FprfAsigID  = EEG.info.DATA(i).fixAsigID(indGrow);
            EEG.info.DATA(i).FprfEventID = EEG.info.DATA(i).fixEventInd(indGrow);
            EEG.info.DATA(i).FprfRank    = EEG.info.DATA(i).fixRank(indGrow);
            
            % Generate First Fixations (FF)
            [~,indUn,~] = unique(EEG.info.DATA(i).FprfAsigW);

            EEG.info.DATA(i).FfAsigW   = EEG.info.DATA(i).FprfAsigW(indUn); 
            EEG.info.DATA(i).FfAsigID  = EEG.info.DATA(i).FprfAsigID(indUn);
            EEG.info.DATA(i).FfEventID = EEG.info.DATA(i).FprfEventID(indUn);
            EEG.info.DATA(i).FfRank    = EEG.info.DATA(i).FprfRank(indUn);
            
            % Generate fixation start time, related to Sentence Start
            fixTimes = [EEG.event(EEG.info.DATA(i).FprfEventID).latency];
            EEG.info.DATA(i).fixTimes  = fixTimes - fixTimes(1);
        end
    end
    function EEG = epochFix(obj)
        EEG = pop_loadset(obj.cfg.inFile);
        EEG = eeg_checkset(EEG);

        for i = 1:length(EEG.info.DATA)
            x = EEG.info.DATA(i);
            EEG.info.DATA(i).SntcLength = repmat(length(x.pos),1,length(x.pos));
        end

        % I only need first fix for each word
        if obj.cfg.firstPass
            fixAsigID    = [EEG.info.DATA.FprfAsigID];
            fixEventInd  = [EEG.info.DATA.FprfEventID];
        else
            fixAsigID    = [EEG.info.DATA.fixAsigID];
            fixEventInd  = [EEG.info.DATA.fixEventInd];
        end
        fixEventRank = [EEG.info.DATA.fixRank];
        fixEventTime = [EEG.info.DATA.fixTimes];

        ID          = [EEG.info.DATA.wrdID];
        pred        = [EEG.info.DATA.pred];
        predPrev    = [EEG.info.DATA.predPrev];
        predNext    = [EEG.info.DATA.predNext];
        sntType     = [EEG.info.DATA.SntType];
        words       = [EEG.info.DATA.palabras];
        SntcLength  = [EEG.info.DATA.SntcLength];
        catgram     = arrayfun(@(x) x.catgram(:)', EEG.info.DATA, 'UniformOutput',0);
        catgram     = [catgram{:}];
        pos         = [EEG.info.DATA.pos];
        lngth       = [EEG.info.DATA.length];
        freq        = [EEG.info.DATA.freq];
        posRelRP    = [EEG.info.DATA.posRelRP];
                
        RP          = ([EEG.info.DATA.posRelRP] == 0) + ([EEG.info.DATA.posRelRP] > 0)*2;
        RP = RP-1;

        [firstFixID, firstFixPos, ~] = unique(fixAsigID);
        firstFixEvent = fixEventInd(firstFixPos);
        fixRank = fixEventRank(firstFixPos); 
        fixTime = fixEventTime(firstFixPos); 
        
        % Remove Baseline and resample to 128Hz
        EEG = pop_epoch(EEG, {}, obj.cfg.timelim, 'eventindices', firstFixEvent);
        EEG = pop_rmbase(EEG, obj.cfg.baseline);
%         EEG = pop_resample(EEG, 128);

        EEG = eeg_checkset(EEG);
        for i = 1:length(EEG.epoch)
            tmp = ID == firstFixID(i);
            EEG.epoch(i).ID      = firstFixID(i);
            EEG.epoch(i).pred    = pred(tmp);
            EEG.epoch(i).predPrev= predPrev(tmp);
            EEG.epoch(i).predNext= predNext(tmp);
            EEG.epoch(i).sntType = sntType(tmp);
            EEG.epoch(i).words   = words(tmp);
            EEG.epoch(i).SntcLength = SntcLength(tmp);
            EEG.epoch(i).catgram = catgram(tmp);
            EEG.epoch(i).pos     = pos(tmp);
            EEG.epoch(i).length  = lngth(tmp);
            EEG.epoch(i).freq    = freq(tmp);
            EEG.epoch(i).RP      = RP(tmp);
            EEG.epoch(i).posRelRP= posRelRP(tmp);
            
            EEG.epoch(i).fixRank = fixRank(i);
            EEG.epoch(i).fixTime = fixTime(i);

            % Find current Fixation Duration
            indF = find(cellfun(@(x) x, EEG.epoch(i).eventlatency) == 0,1);
            EEG.epoch(i).fixDur = EEG.epoch(i).eventduration{indF};

            % Find previous Saccade Duration
            f = regexp({EEG.epoch(i).eventtype{1:indF}}, 'saccade');
            indS = find(cellfun(@(x) ~isempty(x), f),1, 'last');
            if ~isempty(indS)
                EEG.epoch(i).saccDur = EEG.epoch(i).eventduration{indS};
            else
                fprintf('Previous Saccade not found\n')
                EEG.epoch(i).saccDur = NaN;
            end

            % Find previous Fixation Duration
            f = regexp({EEG.epoch(i).eventtype{1:indF-1}}, 'fixation');
            indPF = find(cellfun(@(x) ~isempty(x), f),1, 'last');
            if ~isempty(indPF)
                EEG.epoch(i).prevFixDur = EEG.epoch(i).eventduration{indPF};
            else
                fprintf('Previous Fixation not found\n')
                EEG.epoch(i).prevFixDur = NaN;
            end
        end
    end
    function EEG = generateEvents(obj)
        EEG = pop_loadset(obj.cfg.inFile);
        EEG = eeg_checkset(EEG);
        
        % Inizialize isFixOnWord to 0
        C = num2cell(zeros(1,length([EEG.event.type])));
        [EEG.event.isFixOnWord] = C{:};
        
        % Expand word Properties (generate some new ones)
        for i = 1:length(EEG.info.DATA)
            EEG = addFirstFix(EEG,i);
            x = EEG.info.DATA(i);
            EEG.info.DATA(i).SntcLength = repmat(length(x.pos),1,length(x.pos));
        end
        notFirstFix = ~[EEG.info.DATA.isFirstFix]*1;
        ID          = [EEG.info.DATA.wrdID];
        pred        = zscore([EEG.info.DATA.pred]);
        predPrev    = zscore([EEG.info.DATA.predPrev]);
        predNext    = zscore([EEG.info.DATA.predNext]);
        sntType0    = ~[EEG.info.DATA.SntType]*1;
        sntType1    = [EEG.info.DATA.SntType];
        words       = [EEG.info.DATA.palabras];
        SntcLength  = [EEG.info.DATA.SntcLength];
        catgram     = arrayfun(@(x) x.catgram(:)', EEG.info.DATA, 'UniformOutput',0);
        catgram     = [catgram{:}];
        pos         = zscore([EEG.info.DATA.pos]);
        lngth       = zscore([EEG.info.DATA.length]);
        freq        = zscore([EEG.info.DATA.freq]);
        posRelRP    = [EEG.info.DATA.posRelRP];
        preRP       = ([EEG.info.DATA.posRelRP]< 0)*1;
        postRP      = ([EEG.info.DATA.posRelRP]>=0)*1;
        FprfRank    = ([EEG.info.DATA.FprfRank]);
        
        keyboard
        % pre-RP = -1; Rp = 0; post-RP = 1
        RP          = ([EEG.info.DATA.posRelRP] == 0) + ([EEG.info.DATA.posRelRP] > 0)*2;
        RP = RP-1;

        sujNum = repmat(obj.cfg.iSuj,1,length(RP));

        % Expand fixation asignation
        fixAsigID   = [EEG.info.DATA.fixAsigID];
        fixEventInd = [EEG.info.DATA.fixEventInd];
        
        % Copy information from DATA to event, by Fixation IDs
        for i = 1:length(fixAsigID)
%             fprintf('%d / %d \n',i,length(fixAsigID))
            i_event = fixEventInd(i);
            i_stim  = ID == fixAsigID(i);
            
            EEG.event(i_event).isFixOnWord = 1;
            % Word Properties
            EEG.event(i_event).ID       = ID(i_stim);
            EEG.event(i_event).pred     = pred(i_stim);
            EEG.event(i_event).predPrev = predPrev(i_stim);
            EEG.event(i_event).predNext = predNext(i_stim);
            EEG.event(i_event).sntType0  = sntType0(i_stim);
            EEG.event(i_event).sntType1  = sntType1(i_stim);
            EEG.event(i_event).words    = words(i_stim);
            EEG.event(i_event).catgram  = catgram(i_stim);
            EEG.event(i_event).pos      = pos(i_stim);
            EEG.event(i_event).length   = lngth(i_stim);
            EEG.event(i_event).freq     = freq(i_stim);
            EEG.event(i_event).RP       = RP(i_stim);
            EEG.event(i_event).posRelRP = posRelRP(i_stim);
            EEG.event(i_event).preRP    = preRP(i_stim);
            EEG.event(i_event).postRP   = postRP(i_stim);
            EEG.event(i_event).SntcLength = SntcLength(i_stim);

            EEG.event(i_event).sujNum   = sujNum(i_stim);
            
            EEG.event(i_event).notFirstFix = notFirstFix(i);

            % Find previous Fixation Duration 
            dir = -1; string = 'B_fixation';
            ind = findEvent(EEG, string, i_event, dir);
            EEG.event(i_event).fixDur_prev = EEG.event(ind).duration;

            % Find current Fixation Duration and rank
            ind = i_event;
            EEG.event(i_event).fixDur = EEG.event(ind).duration;
                        
            % Find next Fixation Duration 
            dir = 1; string = 'B_fixation';
            ind = findEvent(EEG, string, i_event, dir);
            EEG.event(i_event).fixDur_next = EEG.event(ind).duration;

            % Find incoming Saccade Duration and Amplitude
            dir = -1; string = 'B_saccade';
            ind = findEvent(EEG, string, i_event, dir);
            EEG.event(i_event).saccDur_in = EEG.event(ind).duration;
            EEG.event(i_event).saccAmp_in = EEG.event(ind).sac_amplitude;
            
            % Find outgoing Saccade Duration and Amplitude
            dir = 1; string = 'B_saccade';
            ind = findEvent(EEG, string, i_event, dir);
            EEG.event(i_event).saccDur_out = EEG.event(ind).duration;
            EEG.event(i_event).saccAmp_out = EEG.event(ind).sac_amplitude;
        end
    end
    function EEG = filterband(obj)
        EEG = pop_loadset(obj.cfg.inFile);
        EEG = eeg_checkset(EEG);
%         EEG = pop_eegfilt(EEG, obj.cfg.locutoff, obj.cfg.hicutoff);
        minphase = true; %filtro causal
        EEG = pop_eegfiltnew(EEG, obj.cfg.locutoff, obj.cfg.hicutoff,[],false,[],[],minphase);
        for ch=1:size(EEG.data,1)
            EEG.data(ch,:,:) = abs(hilbert(EEG.data(ch,:,:)));
        end
    end
end

methods(Static)
    function names = loadSubjects(folder, expe, pre, type)
        % Load filenames (*.type) from folder.
        % expe: [bool] if 1 load suj.$type
        % pre:  [bool] if 1 load suj_pre*.$type
        % type: file extension
        
        names = dir(folder); 
        names = {names(:).name};
        names = names(1,3:end);
        names = names(cellfun(@(x) ~isempty(x),regexp(names, type)));

        if pre == 1 && expe == 1
            return
        elseif pre == 1 
            names = names(cellfun(@(x) ~isempty(x),regexp(names, '_pre')));
        elseif pre == 11
            names = names(cellfun(@(x) ~isempty(x),regexp(names, '_pre1')));
        elseif pre == 12 
            names = names(cellfun(@(x) ~isempty(x),regexp(names, '_pre2')));
        elseif expe == 1 
            names = names(cellfun(@(x) isempty(x),regexp(names, '_pre')));
        end  
    end    
    function badchanslist = badchannels()
        badchanslist = {'ag1',      [26, 34]; ... visual: 34 (Hay un canal jodido que no se si es el B2 u otro), espectro: 26
                'ag1_pre1', [26, 34]; ... visual: 26-34, espectro: 26
                'ag1_pre2', [58, 81]; ... visual: {}, espectro: 81 - 58        
                'at1',      [];...        visual: {}, espectro: {}
                'at1_pre1', [53];...      visual: 53, espectro: 53
                'at1_pre2', [];...        visual: {}, espectro: {}
                'ch1',      [69, 76];...  visual: {}, espectro: 69 - 76
                'ch1_pre1', [69, 76];...  visual: 88, espectro: 88;
                'ch1_pre2', [];...        visual: {}, espectro: {}
                'cm1',      [90];...      visual: 90, espectro: {}
                'cm1_pre1', [];...        visual: {}, espectro: {}
                'cm1_pre2', [];...        visual: {}, espectro: {}
                'cr1',      [104 114 119];... visual: 104 114, espectro: 104 119
                'cr1_pre1', [54  119];... visual: {}, espectro: 54 119
                'cr1_pre2', [115 119];... visual: {}, espectro: 115 119
                'fb1',      [75];...      visual: {}, espectro: 75
                'fb1_pre1', [75];...      visual: {}, espectro: 75
                'fb1_pre2', [75];...      visual: {}, espectro: 75
                'fr1',      [3];...       visual: 3,  espectro: 3
                'fr1_pre1', [];...        visual: {}, espectro: {}
                'fr1_pre2', [3];...       visual: 3,  espectro: 3
                'ga1',      [];...        visual: {}, espectro: {}
                'ga1_pre1', [113];...     visual: {}, espectro: 113
                'ga1_pre2', [];...        visual: {}, espectro: {}
                'gm1',      [];...        visual: {}, espectro: {}
                'gm1_pre1', [];...        visual: {}, espectro: {}
                'gm1_pre2', [];...        visual: {}, espectro: {}
                'ib1',      [92 93 94];...visual: 92 93 94, espectro: {}
                'ib1_pre1', [92 93 94];...visual: 92 93 94, espectro: {}
                'ib1_pre2', [92 93 94 114 116];...visual: 92 93 94, espectro: 114 116
                'lf1',      [53 94 95];...visual: 53 94 95, espectro: {}
                'lf1_pre1', [53 94 95];...visual: 53 94 95, espectro: {}
                'lf1_pre2', [53 94 95];...visual: 53 94 95, espectro: {}
                'lg1',      [52 53 54 64];...     visual: {}, espectro: 52 53 54 64
                'lg1_pre1', [92 93];...   visual: {}, espectro: 92 93
                'lg1_pre2', [54 ];...     visual: {}, espectro: 54
                'mk1',      [12 13 107];...visual: 12 13, espectro: 107
                'mk1_pre1', [];...        visual: {}, espectro: {}
                'mk1_pre2', [26];...      visual: {}, espectro: 26
                'mp1',      [7];...       visual: {}, espectro: 7    
                'mp1_pre1', [7 80];...    visual: 7,  espectro: 80
                'mp1_pre2', [80];...      visual: {}, espectro: 80
                'ng1',      [52 80];...   visual: 80, espectro: 52
                'ng1_pre1', [80];...      visual: {}, espectro: 80
                'ng1_pre2', [80];...      visual: {}, espectro: 80
                'gj1',      [71];...      visual: {}, espectro: 71
                'gj1_pre1', [71];...      visual: {}, espectro: 71
                'gj1_pre2', [71];...      visual: {}, espectro: 71  %OJO volver a ver dsp de interpolar que parece haber otro electrodo mal
                'gw1',      [115 124];... visual: {}, espectro: 115, 124
                'gw1_pre1', [72 82 83 92 94 115 124];...   visual: {}, espectro: 72 82 83 92 94 115 124
                'gw1_pre2', [94 124];...  visual: {}, espectro: 94 124
                'je1',      [103];...     visual: {}, espectro: 103
                'je1_pre1', [55 60 62 101];... visual: 70, espectro: 55 60 62 101
                'je1_pre2', [12 62 71];...visual: {}, espectro: 12 62 71
                'ld1',      [];...        visual: {}, espectro: 
                'ld1_pre1', [24 25 41 42];...              visual: {}, espectro: 24 25 41 42 
                'ld1_pre2', [25 104];...  visual: {}, espectro: 25 104
                'ls1',      [22 80];...   visual: {}, espectro: 22 80
                'ls1_pre1', [];...        visual: {}, espectro: 
                'ls1_pre2', [103 105 106];... visual: {}, espectro: 103 105 106 
                'ml2',      [82 92 93];...    visual: {}, espectro: 82 92 93
                'ml2_pre1', [72 82 92 93];... visual: {}, espectro: 72 82 92 93
                'ml2_pre2', [69 74 93 ];...   visual: {}, espectro: 
                'mo1',      [41 60 70 103];...             visual: {}, espectro: 60 103 
                'mo1_pre1', [59 105 106 118 119];...       visual: {}, espectro: 59 105 106 118 119
                'mo1_pre2', [41 60];...                    visual: {}, espectro: 41 60
                'rl1',      [18 71 81];...             visual: {}, espectro: 18 71 81
                'rl1_pre1', [71 81 119];...            visual: {}, espectro: 71 81 119
                'rl1_pre2', [18 68 70 71 81];...       visual: {}, espectro: 18 68 70 71 81
                'sd1',      [];...                  visual: {}, espectro: 
                'sd1_pre1', [15 16 17 18 31 32];... visual: {}, espectro: 15 16 17 18 31 32
                'sd1_pre2', [41 59];...             visual: 41, espectro: 
                'ci1', [53];...               visual: , espectro: 53 
                'ci1_pre1', [19,53, 110];...  visual: , espectro: 19,53, 110 
                'ci1_pre2', [53,104];...      visual: , espectro: 53,104 
                'dc1', [58,69];...                 visual: , espectro: 58,69
                'dc1_pre1', [58,69];...            visual: , espectro: 58,69
                'dc1_pre2', [14,72];...            visual: , espectro: 14,72
                'df1', [12, 71, 87];...               visual: , espectro: 12, 71, 87
                'df1_pre1', [10,11,12,15];...         visual: , espectro: 10,11,12,15
                'df1_pre2', [12,71,72];...            visual: , espectro: 12,71,72
                'gc1', [];...                  visual: , espectro: 
                'gc1_pre1', [68,108,109];...   visual: , espectro: 68,108,109
                'gc1_pre2', [67,68];...        visual: , espectro: 67,68
                'gl1', [38];...                     visual: , espectro: 38
                'gl1_pre1', [38,44,118];...         visual: , espectro: 38,44,118
                'gl1_pre2', [38,101];...            visual: , espectro: 38,101
                'il1', [9,55];...                       visual: , espectro: 9, 55
                'il1_pre1', [12, 55, 63];...            visual: , espectro: 12, 55, 63
                'il1_pre2', [54, 55, 63];...            visual: , espectro: 54, 55, 63
                'jc1', [58, 95];...                visual: , espectro: 58,95
                'jc1_pre1', [27];...               visual: , espectro: 27
                'jc1_pre2', [57,58];...            visual: , espectro: 57,58
                'va1', [];...                 visual: , espectro: 
                'va1_pre1', [72,73,74,80,81];...            visual: , espectro: 72,73,74,80,81
                'va1_pre2', [];...            visual: , espectro: 
                };
    end
    function marks = syncMarks(sujName)
        marks = [254 255];
        
        if ismember(sujName, {'ag1_pre1','cm1_pre1','cm1_pre2','cr1_pre1','gm1_pre1','lf1_pre1','mp1_pre1','ng1_pre1','gj1_pre1','gw1_pre1', 'je1_pre1', 'ld1_pre1', 'ml2_pre1', 'mo1_pre1', 'sd1_pre1', 'gc1_pre1', 'gl1_pre1', 'jc1_pre1', 'va1_pre1'})
            marks = [100 255];
        elseif ismember(sujName, {'ls1_pre1','rl1_pre2','il1_pre1'})
            marks = [101 255];
        elseif strcmpi(sujName, 'gm1')
            marks = [112 255];
        elseif strcmpi(sujName, 'mo1')
            marks = [12 255];
        elseif strcmpi(sujName, 'gl1')
            marks = [3 255];
        end
    end
    function [EEG, ET] = renameEEGevents(EEG, ET)
       if EEG.event(1).type == 255
            EEG.event(1).type = 254;
       end
       
       % Delete nonimportant events
       ET.event = ET.event(ET.event(:,2) > 99,:);
       ind = cellfun(@(x) ~ismember(x, [100,101,254,255]), {EEG.event.type});
       EEG.event(ind) = [];
    end
    function EEG = generateBestEyeChan(EEG)
        sujName = EEG.filename(1:end-4);
        if strcmpi(sujName, 'gm1')
            keyboard
        elseif strcmpi(sujName, 'mo1')
            keyboard
        end
        
        % split in blocks
        % for some reason there are more 255 marks than there should be
        % i look for all the 254 and 255 marks
        % during rests there are no EEG marks, only ET marks
        % I want to delete rests moments
        ET = load(['0-load_ET/' sujName '.mat']);  
        cals  = findCalibratios(ET.messages);
        
        events = {EEG.event.type};
        events = arrayfun(@(x) str2double(x),events); 
        times  = [EEG.event.latency];        
%         scatter(times,events)

        if ~isempty(strfind(sujName, '_pre'))
            tStarts = [1]; % from the first sample 
            tEnds   = [length(EEG.data)]; % to the last one
        else
            cals = [cals(1) cals]; % repeate first cal for training (share cal with first block)

            % Find blocks starts
            tStarts = times(events == 254);        
            tStarts(1) = 1; % from the first sample 
            
            % Define blocks ends
            tEnds = [(tStarts(2:end)-1) length(EEG.data)];
        end
        
        ocularChans = [];
        for i=1:length(tStarts)
            if length(EEG.chanlocs) == 148; % monocular
                eye = cals(i).bestCal;
                chans = find(cellfun(@(x) ~isempty(x), strfind({EEG.chanlocs.labels},'-GAZE-X')));
            else
                eye = cals(i).bestCal;
                chans = find(strcmpi([eye '-GAZE-X'], {EEG.chanlocs.labels}));
            end
            
            chans = [chans:chans+2];

            t = tStarts(i):tEnds(i);
            ocularChans = [ocularChans EEG.data(chans,t)];
            
            inds = find(times >= tStarts(i) & times < tEnds(i));
            
%             arrayfun(@(x) strrep(x.type, [eye '_'], ['B_']), EEG.event(inds))
            
            for j=1:length(inds)
                e =  EEG.event(inds(j)).type;           
                EEG.event(inds(j)).type = strrep(e, [eye '_'], ['B_']);
            end
        end

        EEG.chanlocs = [EEG.chanlocs; EEG.chanlocs(chans)];
        EEG.chanlocs(end-2).labels = 'BEST-GAZE-X';
        EEG.chanlocs(end-1).labels = 'BEST-GAZE-Y';
        EEG.chanlocs(end-0).labels = 'BEST-AREA';
        EEG.nbchan = length(EEG.chanlocs);

        EEG.data = [EEG.data; ocularChans]; 
        EEG = eeg_checkset(EEG);

        % remove L and R channels
        indDlt1 = find(~cellfun(@isempty,regexp({EEG.chanlocs.labels}, '[RL]-')));
        indDlt2 = find(~cellfun(@isempty,regexp({EEG.chanlocs.labels}, 'TIME')));
        indDlt3 = find(~cellfun(@isempty,regexp({EEG.chanlocs.labels}, 'EXG[5-8]')));
        indDlt  = [indDlt1 indDlt2 indDlt3];
        EEG = pop_select(EEG, 'nochannel', indDlt);
        EEG = eeg_checkset( EEG );
    end
    function EEG = prepareDataForICA(inFile, icaFilterEdge,conf)

        % Load EEG        
        EEG = pop_loadset(inFile);
        EEG = eeg_checkset(EEG);  
        % Remove unuse ET channels
        indDlt = find(~cellfun(@isempty,regexp({EEG.chanlocs.labels}, '[RL]-')));
        EEG = pop_select(EEG, 'nochannel', indDlt);
        EEG = eeg_checkset( EEG );
        % Filter
        tmp = EEG.data(129:end,:);
        EEG.data = EEG.data(1:128,:);
        EEG = pop_eegfiltnew(EEG, icaFilterEdge,[]);      
        EEG.data = [EEG.data; tmp];
        EEG = eeg_checkset( EEG );
        % Only keep sentences
        events = {EEG.event.type};
        events = arrayfun(@(x) str2double(x),events); 
        times  = [EEG.event.latency];        
        senteceStart = times(events==220);
        senteceEnds  = times(events==221);
        if strcmpi(EEG.filename, 'gl1.set')
            senteceEnds  = senteceEnds(2:end);
        elseif strcmpi(EEG.filename, 'gm1.set')
            keyboard % senteceEnds  = senteceEnds(2:end);
        elseif strcmpi(EEG.filename, 'mo1.set')
            keyboard % senteceEnds  = senteceEnds(2:end);
        end
        
        ind  = [senteceStart'  senteceEnds'];
        EEG = pop_select(EEG, 'nopoint', ind); 
        EEG = eeg_checkset( EEG );
        
        % Load pre1
        pre1 = pop_loadset([inFile(1:end-4) '_pre1.set']);
        % Remove unuse ET channels
        indDlt = find(~cellfun(@isempty,regexp({pre1.chanlocs.labels}, '[RL]-')));
        pre1 = pop_select(pre1, 'nochannel', indDlt);
        % Filter
        tmp = pre1.data(129:end,:);
        pre1.data = pre1.data(1:128,:);
        pre1 = pop_eegfiltnew(pre1, icaFilterEdge,[]);      
        pre1.data = [pre1.data; tmp];
        pre1 = eeg_checkset(pre1);
        % keep only excersices
        events = {pre1.event.type};
        events = arrayfun(@(x) str2double(x),events); 
        times  = [pre1.event.latency];        
        starts = times(events==100);
        ends  = times(events==101);
        if strcmpi(pre1.filename, 'il1_pre1.set')
            ends  = ends(2:end);
        end
        ind  = [starts'  ends'];
        pre1 = pop_select(pre1, 'nopoint', ind); 
        pre1 = eeg_checkset( pre1 );
                
        % Load pre2
        pre2 = pop_loadset([inFile(1:end-4) '_pre2.set']);
        % Remove unuse ET channels
        indDlt = find(~cellfun(@isempty,regexp({pre2.chanlocs.labels}, '[RL]-')));
        pre2 = pop_select(pre2, 'nochannel', indDlt);
        % Filter
        tmp = pre2.data(129:end,:);
        pre2.data = pre2.data(1:128,:);
        pre2 = pop_eegfiltnew(pre2, icaFilterEdge,[]);      
        pre2.data = [pre2.data; tmp];
        pre2 = eeg_checkset(pre2);
        % keep only excersices
        events = {pre2.event.type};
        events = arrayfun(@(x) str2double(x),events); 
        times  = [pre2.event.latency];        
        starts = times(events==100);
        ends  = times(events==101);
        ind  = [starts'  ends'];
        pre2 = pop_select(pre2, 'nopoint', ind); 
        pre2 = eeg_checkset( pre2 );
        
        % mergre Pre1 and Pre2
        PRE = pop_mergeset(pre1, pre2);
        PRE = eeg_checkset( PRE );
        
        % Merge PRE and EEG
        EEG = pop_mergeset(PRE, EEG);
        EEG = eeg_checkset(EEG);
    end
    function EEG = preparePre(EEG, inFile)
        
        % Remove unuse ET channels
        indDlt = find(~cellfun(@isempty,regexp({EEG.chanlocs.labels}, '[RL]-')));
        PRE1 = pop_select(EEG, 'nochannel', indDlt);
        PRE1 = eeg_checkset(PRE1);      
        
        fname = [inFile(1:end-5) '2.set'];
        
        PRE2 = pop_loadset(fname);
        indDlt = find(~cellfun(@isempty,regexp({PRE2.chanlocs.labels}, '[RL]-')));
        PRE2 = pop_select(PRE2, 'nochannel', indDlt);
        PRE2 = eeg_checkset(PRE2);     

        % mergre Pre1 and Pre2
        EEG = pop_mergeset(PRE1, PRE2);
        EEG = eeg_checkset(EEG);
        
    end
    function EEG = addChanInfo(EEG)
        
        EEG = pop_select( EEG,'nochannel',129:EEG.nbchan);
        EEG = eeg_checkset( EEG );

        load coords_LNI_128_toreplaceinEEG
        EEG.chanlocs    = CHANS.chanlocs;
        EEG.urchanlocs  = CHANS.urchanlocs;
        EEG.chaninfo    = CHANS.chaninfo;
        EEG = eeg_checkset( EEG );

    end
    function [EEG, new_order] = shuffleEvents(EEG, fields, within)
        new_order = [];
        if within
            [C,ia,ic]=unique([EEG.event.sujNum]);
            for i = 1:length(C)
                start = find(ic==i,1,'first');
                stop  = find(ic==i,1,'last');
                new_order = [new_order randperm(stop-start+1)+start-1];
            end
        else
            new_order = randperm(length(EEG.event));
        end
        
        for i=1:length(fields)
            c = fields{i};
            [EEG.event(:).(c)] = EEG.event(new_order).(c);
        end

    end
    
end%methods
end% classdef


function cals  = findCalibratios(messages)
    fprintf('Finding calibrations marks\n')

    C = messages;
    [stLinesNum, ~]           = findText(C, {'CALIBRATION'}, 0);
    [modeLinesNum, modeLines] = findText(C, {'!MODE'}, 0);

    headerLim  = [];
    % For each block find the last calibration
    for i = 1:length(modeLinesNum)
        mode     = strsplit(strtrim(modeLines{i}), ' '); 
        mode     = mode(length(mode));
        previous = find((stLinesNum  - modeLinesNum(i)) < 0);

        if strcmpi(mode, 'L') || strcmpi(mode, 'R')
            % If MONO, I only care for the last calibration
            calStart = stLinesNum(length(previous));
        elseif strcmpi(mode, 'LR')
            % If BINO, I need the last 2 calibrations (L and R)
            calStart = stLinesNum(length(previous)-1);
        end

        headerLim  = [headerLim  ; [calStart modeLinesNum(i)]];
    end

    % if there was a fail start there is 2 calibrations for 1 block
    diffStarts = diff(headerLim(:,1));
    index      = find(diffStarts == 0);
    headerLim(index,:)  = [];

    for i = 1:size(headerLim,1)            
        tCal = C(headerLim(i,1) : headerLim(i,2));
        cals(i) = findBestEye(tCal);
    end
end
function [line, lines] = findText(C, texts, first)
    % find one of the texts in file from line 1 to line maxlines
    % output: # of lines and text of line

    line  = [];
    lines = [];
    for i = 1:length(texts)% find each of str in texts
        if first == 1
            counter = 0;
%             keyboard
            while isempty(line) && counter < length(C)
                counter = counter+1;
                tline = C{counter};    
                if ~isempty(strfind(tline,texts{i})) 
                    line = [line counter];
                    lines = [lines tline];
                    return 
                end
                %    disp([num2str(contador) tline])
            end
        else
            a = cellfun(@(x) ~isempty(strfind(x,texts{i})), C);
            index = find(a);
            line  = [line index];
            lines = [lines C(line)];
        end                

    end
%     if isempty(line)
%         tline = [];
%     end
end
function all   = findBestEye(cal)

    fprintf('Finding best eye:\n')
    [~, lines] = findText(cal, {'!CAL VALIDATION'},0);
    
    splited = strsplit(lines{end});
    L = any(findText(splited, {'LEFT'},0));
    R = any(findText(splited, {'RIGHT'},0));
    
    if xor(L, R) 
        eyes = {'L' 'R'};
        all.eye     = eyes{[L R]};
        all.bestCal = all.eye;

        ind = find(cellfun(@(x) ~isempty(x), (strfind(splited,'ERROR'))));
        all.calErr  = str2num(splited{ind+1});
        
        fprintf('\t%s: %.2f\n',  all.eye,...
                                 all.calErr)
        fprintf('\tBest eye: %s\n', all.bestCal)

    elseif L && R
        all.eye = 'BOTH';
        
        indR = find(cellfun(@(x) ~isempty(x), (strfind(splited,'ERROR'))));
        abortedR = any(cellfun(@(x) ~isempty(x), (strfind(splited,'ABORTED'))));
        if abortedR; calR = inf; else; calR = str2num(splited{indR+1}); end

        indL = find(cellfun(@(x) ~isempty(x), (strfind(splited,'ERROR'))));
        abortedL = any(cellfun(@(x) ~isempty(x), (strfind(splited,'ABORTED'))));
        if abortedL; calL = inf; else; calL = str2num(splited{indL+1}); end

        if calR < calL
            all.bestCal = 'R';
            all.calErr  = calR;
        else
            all.bestCal = 'L';
            all.calErr  = calL;
        end


        fprintf('\tRIGHT: %.2f  LEFT: %.2f\n',  calR, calL)
        fprintf('\tBest eye: %s\n', all.bestCal)

    else 
        keyboard
        fprintf('Unknwon eye\n')
    end
end
function ppc   = calculatePPC(DATA)
    % length of sentence image / length (in chars) of sentence
    % mean of all the sentences
    tmp = [];
    for i = 1:length(DATA)
        tmp = [tmp size(DATA(1).imagen,2)/length(DATA(1).oracion)];
    end

    ppc = mean(tmp); % Pixel per character
end
function DATA  = addFieldsEstim(DATA,ESTIM,campos)
    try
        fprintf('\tAdd fields from ESTIM to DATA\n')
        for i=1:length(DATA)
            for indcampo=1:length(campos)
                DATA(i).(campos{indcampo})=ESTIM(DATA(i).ind+1).(campos{indcampo});
            end
        end
    catch ME
        ME
        keyboard
    end
end    
function DATA  = addSpaces(DATA)
    % Pixel for character
    ppc = calculatePPC(DATA);
    for trial=1:length(DATA)
        spaces = strfind(DATA(trial).oracion,' ');
        spaces = [0 spaces length(DATA(trial).oracion)+1];        
        if strcmp(DATA(trial).oracion(end),' ')
            spaces(end) = [];
        end
        str = DATA(trial).posoracion(1);
        DATA(trial).spaces = str +(spaces-.5)*ppc;
    end

%             Max Vertical Pos for fixations
%             mvp = init.height/2 + init.height*.2;
%             all = ET.parsed(block);
%             all = interpolateBlinksSaccades(all);            
%             fprintf('\tBest eye: %s\n\n', all.bestCal)                      
%             DATA = addFix(all,DATA,mvp, ppc);
%             DATA = assignFixToWords(DATA, init);
%             DATA = assignFixToChar(DATA, ppc);
%             DATA = assignFixToFirstChar(all, DATA);
%             DATA = addFixCross(all, DATA, init);
%             DATA = addPupil(all, DATA, SUJ(iSuj));
%             DATA = addMsgsToTrial(all, DATA);    
%             DATA = addRT(all, DATA);
%             DATA = addSpaces(DATA, ppc);

end
function DATA  = addLatency(DATA, EEG)
    starts = arrayfun(@(x) ~isempty(x{1}),regexp({EEG.event.type}, '220')); 
    starts = [EEG.event(starts).latency]/EEG.srate;

    ends = arrayfun(@(x) ~isempty(x{1}),regexp({EEG.event.type}, '221')); 
    ends = [EEG.event(ends).latency]/EEG.srate;

    % If recording stated late and the first mark is an end mark
    if length(starts)+1 == length(ends)
        ends = ends(2:end);
    end

    inds  = [starts(4:end)'  ends(4:end)']; % Delete 3 sentences of practice 
    for iData = 1:length(DATA)
        DATA(iData).latency = inds(iData,:);
    end
end
function ind = findEvent(EEG, string, from, dir)
    j = 0;
    found = 0;
    while ~found
        j= j+1;
        found = strcmpi(EEG.event(from+j*dir).type, string);
    end
    ind = from + j*dir;
end
function EEG = addFirstFix(EEG,i)
    f = EEG.info.DATA(i).fixAsigWrd;
    EEG.info.DATA(i).isFirstFix  = zeros(1,length(f));
    for j = 1:length(f)
        EEG.info.DATA(i).isFirstFix(j) = find(f(j) == f,1,'first') == j;
    end
end