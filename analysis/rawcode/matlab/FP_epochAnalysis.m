classdef FP_epochAnalysis
    % FP_epochAnalysis
    %
    % Class containning all functions needed for creating epochs acording
    % to differenct conditions.
    %           
    %
    % properties:
    %           cfg (configuration struct)
    %           Example: 
    %                   inFolder                =
    %                   session_path.data_analysis;fhilb
    %                   fp.cfg.subjname         =  session_path.subjname;
    %                   fp.cfg.sessionfilenames = session_path.sessionfilenames;
    %                   fp.cfg.matdir           = matdir;
    %
    %           Methods:
    %                   fixationsAnalysis(obj)                   
    %                   fixationsAnalysisOne(obj)
    %
    %           Static Methods:
    %                   loadSubjects(folder, expe, pre, type, keyw)      
    %                   fun_detect_succesive_fixs(eyetrial)
    %                   fun_extract_attributes_target(SingleTrial)
    %                   fun_extract_dist_attributes(SingleTrial, rel_position)
    %                   fun_parse_eyedata_EEG_beta(EEG,ExpTrials)
    %                   fun_makeStructToERP(ERPstruct,condVec,iSuj,folder,subjName)
    %                   indexCondEpochs(EEG, fixsSuj,minfixdur)
    %                   indexCondEpochs2(EEG, fixsSuj,minfixdur)
    %                   indexCondEpochs3(EEG, fixsSuj,minfixdur)
    %                   indexCondEpochs4(EEG, fixsSuj,minfixdur)
    %                   indexCondEpochs5(EEG, fixsSuj,minfixdur)
    %                   ept_TFCE(DataFile1, DataFile2, ElecFile, varargin)
    %                   fun_conditionHistogramInd(EEG,fixsSuj,minfixdur)
    %                   renameEventsUnfold(EEG, fixsSuj,minfixdur,istarget) 
    %   The rest of the functions difinitions are functions used by
    %   ept_TFCE.
    %
    % 2021-21-01 Damian Care
    
    properties
        cfg 
    end
    methods
        function obj                                    = FP_epochAnalysis()
%             obj.cfg.inFolder = [];
%             obj.cfg.subjname = [];
%             obj.cfg.sessionfilenames = [];
%             obj.cfg.matdir           =  [];

        end
        function [EEG,fixs,trials,ExpTrials]            = fixationsAnalysisOne(obj,minDurFix,freqFilteredEEG,fix2stim)
            %             input:
            %                 obj: instantiated class (obj.cfg.infile)
            %                                         (obj.cfg.sessiongilenames)
            %                                         (obj.cfg.subjname)
            %                 minDurFix: minimun fix duration % mdf
            %                 freqFilteredEEG: 1 if EEG is filtered by
            %                 freq-band, 0 otherwise % mdf 14/6/21
            %                 fix2stim: 1 for filtering only fixs to stims,
            %                 0 stim + out of stim % mdf 17/6/21
            %
            %             output:
            %                    EEG: EEGLAB struct.
            %
            %                    fixs: struct with all relevant data about
            %                          each fixation run for subject specified on
            %                          cfg.
            %
            %                    trials: struct with all relevant data about each trial run
            %                            for subject specified on cfg.
            %
            %                    ExpTrials: experimetal data.
            %
            %             children: 
            % 2021-21-01 Damian Care
            % 2021-7-06 modified by Maria da Fonseca                           
            %
            inFile              = obj.cfg.inFile
            sessionfilenames    = obj.cfg.sessionfilenames;
            subjname            = obj.cfg.subjname;
            matdir              = obj.cfg.matdir;
            
            resp_cell_VS        = {'sin resp','A','P','X'};
            resp_cell_EX_gender = {'sin resp','male','female','X'};
            resp_cell_EX_obj    = {'sin resp','small','big','X'};
            fixs = [];
            trials = [];
            tmpfix.refix     =  NaN;
            name = printLoadingFile([inFile]);
                  
            %filename = [ inFolder names{iSuj}];
            if freqFilteredEEG
                suj = name(1:3);
            else
                suj = name(7:end-13);
            end
            ind_matfile = find(ismember(subjname,suj));       

            fprintf('Subject %s\n',suj)
            sessionfilename    = sessionfilenames{ind_matfile};

            matNames = obj.loadSubjects(matdir, 1, 1, 'mat');
            if sum(ismember(matNames,[sessionfilename '.mat']))==0
                fprintf('mat file for subject %s was not found\n', suj)
            else
                fprintf('%s: Loading the new raw MAT file %s from %s folder\n',suj,sessionfilename,matdir)
                load(fullfile(matdir,[sessionfilename '.mat']));
                    if suj == 'E03' %repair
                         ExpTrials = ExpTrials(11:end);
                    elseif suj == 'J04'
                        ExpTrials = ExpTrials([1:71 73:95 97:191 193:215]);
                    elseif suj == 'J09'
                        ExpTrials = ExpTrials([1:23 49:215]);
                    elseif suj == 'J10'
                        ExpTrials = ExpTrials([1:95 98:215]); %10648 urevent from missing trial so im changin ([1:95 97:215]) cause first trial secondSEgment got lost 
                    elseif suj == 'J08'
                        ExpTrials = ExpTrials([1:23 25:167]);
                    elseif suj == 'J07'
                        ExpTrials = ExpTrials([1:23 25:215]);
                    elseif suj == 'J05'
                        ExpTrials = ExpTrials([ 1:47 49:95 97:143  145:215]);
                    end
            end
            EEG = pop_loadset(inFile);
            % Keep 'saccade' and 'fixation' events only
            event = EEG.event;
            %me quedo con la deteccion de Engbert
            event(ismember({event.type},{'R_blink','R_saccade','bad_ET','R_fixation','R_saccade'})) = []; 
            eyemap_events = {'211' '212' '221' '222' '231' '232' '241' '242' '251' '252' '290' '300'};
            event(ismember({event.type},eyemap_events)) = [];
            other_events = {'1' '11' '12' '15' '16' '17' '100' '101' '150' '160' '170' '180' '200'};
            event(ismember({event.type},other_events)) = [];
            %unique({event.type})
            eyedata = obj.fun_parse_eyedata_EEG_beta(EEG,ExpTrials);
            %keyboard
            for tr = 1:length(eyedata) %trials
                C = NaN;%agregado para intentar arreglar J03 etc

                if eyedata(tr).Nfix>0
                    %try
                    [cat,atr,resp]= obj.fun_extract_attributes_target(ExpTrials(tr));
                    %catch e
                    %    fprintf(e.identifier)
                        %keyboard
                    %end
                    if strcmp(ExpTrials(tr).info.trialtype,cfg.tasks{1});%(VS)
                           C = strcmp(resp_cell_VS{ExpTrials(tr).resp.key+1} , resp);
                    elseif strcmp(ExpTrials(tr).info.trialtype,cfg.tasks{2})%(EX)
                        if strcmp(cat,'objects');
                           C = strcmp(resp_cell_EX_obj{ExpTrials(tr).resp.key+1} , atr);
                        elseif strcmp(cat,'faces');
                           C = strcmp(resp_cell_EX_gender{ExpTrials(tr).resp.key+1} , atr);
                        end
                    end

                    % Struct fixs (all subjects)
                    [fix2keep ind_suc mask] = obj.select_multifix(eyedata(tr),minDurFix);
                                                   %crear ayedata.fixs_added y ayadata.rel pos y abspos
                                                   
                    [non_succesive_mask]  = obj.select_non_succesive_fixs(eyedata(tr),minDurFix);      
                    if sum(mask'.*non_succesive_mask)
                        fprintf('the multifix set and the non-multifix set are overlapped');
                    end         
                    
                    % firstFix2target: n value to the first fix to target
                    % that last more than 200 ms
                    firstFix2target = NaN;
                    for n = 1:eyedata(tr).Nfix  %fixations for tr and suj
                        dur     = eyedata(tr).fixs(n,2) - eyedata(tr).fixs(n,1); % mdf 3/6/21 fix duration from E
                        if eyedata(tr).fix2target(n) & dur > 0.2 & isnan(firstFix2target)
                            firstFix2target = n;
                        end
                    end
                    
                    rank = 1;
                    fixbypass = 1;
                    for n = 1:eyedata(tr).Nfix  %fixations for tr and suj
                        tmpdur     = eyedata(tr).fixs(n,2) - eyedata(tr).fixs(n,1); % mdf 3/6/21 fix duration from E
                            
                        if mask(n)
                            fixbypass = 0;
                            multifix = 1;
                        elseif non_succesive_mask(n)
                            fixbypass = 0;    
                            multifix = 0;
                        else
                            fixbypass = 1;    
                            multifix = 1;                    
                        end 
                        
                        % bypass fixs out of stims or not % mdf 17/6/21
                        fixbypass2 = 1;
                        fixbypass3 = 1;
                        if fix2stim & eyedata(tr).stposabs(n)~=0
                            fixbypass2 = 0;
                        elseif fix2stim & eyedata(tr).stposabs(n)==0
                            fixbypass2 = 1;
                        else
                            fixbypass3 = 0;
                        end 
                        
                        if ~fixbypass & (~fixbypass2 | ~fixbypass3)
                            if ~isnan(firstFix2target)
                                if n<firstFix2target
                                    tmpfix.pretarget = 1;
                                else
                                    tmpfix.pretarget = 0;
                                end
                            else
                                tmpfix.pretarget = 1;
                            end
                            tmpfix.multifix = multifix;
                            tmpfix.urevent = eyedata(tr).fixs(n,7);   
                            tmpfix.tonset  = eyedata(tr).fixs(n,1) - eyedata(tr).t_bgn_tr;   
                            tmpfix.dur_added = tmpdur; 
                            tmpfix.rank = rank;
                            if rank<5
                                tmpfix.first_rank=1;
                            else
                                tmpfix.first_rank=0;
                            end
                            rank = rank+1;
                            tmpfix.stposrel= eyedata(tr).stposrel(n);
                            tmpfix.stposabs  = eyedata(tr).stposabs(n);
                            tmpfix.n       = n;
                            tmpfix.x       = eyedata(tr).fixs(n,4); 
                            tmpfix.y       = eyedata(tr).fixs(n,5);
                            tmpfix.istarget= eyedata(tr).fix2target(n); 

                            if tmpfix.stposrel~=ExpTrials(tr).info.indtarget & tmpfix.stposabs~=0 ...
                                    & (sum(ismember(ExpTrials(tr).info.indpos,tmpfix.stposabs))>0)
                                
                                [cat_dist,atr_dist] = obj.fun_extract_dist_attributes(ExpTrials(tr), tmpfix.stposabs);
                                                        %esta funcion deberia arrojar cat y atr para distractor

                                tmpfix.isdistractor = 1;
                                tmpfix.distractor_cat          = cat_dist;
                                tmpfix.ditractor_atr           = atr_dist;

                            else
                                tmpfix.isdistractor = 0;
                                tmpfix.distractor_cat          = NaN;
                                tmpfix.ditractor_atr           = NaN;
%                                 if tmpfix.istarget==0
%                                     keyboard
%                                 end
                            end
                            tmpfix.trial_targetseen     = any(eyedata(tr).fix2target);
                            tmpfix.trial_correct        = C;                      
                            tmpfix.trial_cat            = cat;
                            tmpfix.trial_resp           = resp;
                            tmpfix.trial_atr            = atr;
                            tmpfix.trial_type           = ExpTrials(tr).info.trialtype;
                            tmpfix.trial_cattype         = ExpTrials(tr).info.cattype;
                            tmpfix.trial_number         = tr;
                            tmpfix.trial_subject        = suj;
                            tmpfix.trial_bgn_latency    = eyedata(tr).t_bgn_tr;
                            %change

                            fixs = [fixs, tmpfix];
                            
                        end

                    end

                    % Struct trials (all subjects)
                    if any(eyedata(tr).fix2target)
                        tmptrial.nfix2target    = find(eyedata(tr).fix2target,1,'first');
                        %analizar criterio (n al target tal que fixdur > a 
                        %el promedio de fixdur para distractores) )
                    else
                        tmptrial.nfix2target    = NaN;
                    end
                    tmptrial.Nfix           = eyedata(tr).Nfix;
                    tmptrial.targetseen     = any(eyedata(tr).fix2target);
                    tmptrial.correct        = C;
                    tmptrial.cat            = cat;
                    tmptrial.resp           = resp;
                    tmptrial.atr            = atr;
                    tmptrial.type           = ExpTrials(tr).info.trialtype;
                    tmptrial.number         = tr;
                    tmptrial.subject        = suj;
                    tmptrial.valid          = 1;

                    %CREO vector con las refijaciones para definir tmpfix.refix

                    ind         =  find([fixs.trial_number]==tr ...
                        & arrayfun( @(x)strcmp(x,suj),{fixs.trial_subject}));
                    tmp_pos_rel = [fixs(ind).stposrel];
                    tmp_refix   = zeros(size(tmp_pos_rel));

                    for i = 1:15
                        ind_refix = find(tmp_pos_rel==i);
                        if length(ind_refix)>1
                            tmp_refix(ind_refix) =  1:length(ind_refix);
                        end
                    end
                    ref              = num2cell([tmp_refix]);
                   [fixs(ind).refix] = deal(ref{:});
                else
                    tmptrial.valid   = 0;
                end
               
                trials                  = [trials, tmptrial];

            end 
        end
        function [EEG,fixs,trials,ExpTrials]            = fixationsAnalysisAll(obj) % mdf 30/4/21to unpack multifixs
            %             input:
            %                 obj: instantiated class (obj.cfg.infile)
            %                                         (obj.cfg.sessiongilenames)
            %                                         (obj.cfg.subjname)
            %
            %             output:
            %                    EEG: EEGLAB struct.
            %
            %                    fixs: struct with all relevant data about
            %                          each fixation run for subject specified on
            %                          cfg.
            %
            %                    trials: struct with all relevant data about each trial run
            %                            for subject specified on cfg.
            %
            %                    ExpTrials: experimetal data.
            %
            %             children: 
            % 2021-21-01 Damian Care
            %                              
            %
            inFile              = obj.cfg.inFile
            sessionfilenames    = obj.cfg.sessionfilenames;
            subjname            = obj.cfg.subjname;
            matdir              = obj.cfg.matdir;
            resp_cell_VS        = {'sin resp','A','P','X'};
            resp_cell_EX_gender = {'sin resp','male','female','X'};
            resp_cell_EX_obj    = {'sin resp','small','big','X'};
            fixs = [];
            trials = [];
            tmpfix.refix     =  NaN;
            name = printLoadingFile([inFile]);
                  
            %filename = [ inFolder names{iSuj}];
            suj = name(7:end-13);
            ind_matfile = find(ismember(subjname,suj));       
            %  keyboard
            fprintf('Subject %s\n',suj)
            sessionfilename    = sessionfilenames{ind_matfile};

            matNames = obj.loadSubjects(matdir, 1, 1, 'mat');
            if sum(ismember(matNames,[sessionfilename '.mat']))==0
                fprintf('mat file for subject %s was not found\n', suj)
            else
                fprintf('%s: Loading the new raw MAT file %s from %s folder\n',suj,sessionfilename,matdir)
                load(fullfile(matdir,[sessionfilename '.mat']));
                    if suj == 'E03' %repair
                         ExpTrials = ExpTrials(11:end);
                    elseif suj == 'J04'
                        ExpTrials = ExpTrials([1:71 73:95 97:191 193:215]);
                    elseif suj == 'J09'
                        ExpTrials = ExpTrials([1:23 49:215]);
                    elseif suj == 'J10'
                        ExpTrials = ExpTrials([1:95 98:215]); %10648 urevent from missing trial so im changin ([1:95 97:215]) cause first trial secondSEgment got lost 
                    elseif suj == 'J08'
                        ExpTrials = ExpTrials([1:23 25:167]);
                    elseif suj == 'J07'
                        ExpTrials = ExpTrials([1:23 25:215]);
                    elseif suj == 'J05'
                        ExpTrials = ExpTrials([ 1:47 49:95 97:143  145:215]);
                    end
            end
            EEG = pop_loadset(inFile);
            % Keep 'saccade' and 'fixation' events only
            event = EEG.event;
            %me quedo con la deteccion de Engbert
            event(ismember({event.type},{'R_blink','R_saccade','bad_ET','R_fixation','R_saccade'})) = []; 
            eyemap_events = {'211' '212' '221' '222' '231' '232' '241' '242' '251' '252' '290' '300'};
            event(ismember({event.type},eyemap_events)) = [];
            other_events = {'1' '11' '12' '15' '16' '17' '100' '101' '150' '160' '170' '180' '200'};
            event(ismember({event.type},other_events)) = [];
            %unique({event.type})
            eyedata = obj.fun_parse_eyedata_EEG_beta(EEG,ExpTrials);
            %keyboard
            for tr = 1:length(eyedata) %trials
                C = NaN;%agregado para intentar arreglar J03 etc

                if eyedata(tr).Nfix>0
                    %try
                    [cat,atr,resp]= obj.fun_extract_attributes_target(ExpTrials(tr));
                    %catch e
                    %    fprintf(e.identifier)
                        %keyboard
                    %end
                    if strcmp(ExpTrials(tr).info.trialtype,cfg.tasks{1});%(VS)
                           C = strcmp(resp_cell_VS{ExpTrials(tr).resp.key+1} , resp);
                    elseif strcmp(ExpTrials(tr).info.trialtype,cfg.tasks{2})%(EX)
                        if strcmp(cat,'objects');
                           C = strcmp(resp_cell_EX_obj{ExpTrials(tr).resp.key+1} , atr);
                        elseif strcmp(cat,'faces');
                           C = strcmp(resp_cell_EX_gender{ExpTrials(tr).resp.key+1} , atr);
                        end
                    end

                    % Struct fixs (all subjects)
                    % k = 1;
                    k = 0;
                    tmp = 0;
                    tmpdur=[];
                    [vec_succesiv ind_suc mask] = obj.fun_detect_succesive_fixs(eyedata(tr));
                                                   %crear ayedata.fixs_added y ayadata.rel pos y abspos
                    % keyboard;
                    multi_n = 1;
                    for n = 1:eyedata(tr).Nfix  %fixations for tr and suj
                        tmpdur     = eyedata(tr).fixs(n,2) - eyedata(tr).fixs(n,1); 

                        if n<eyedata(tr).Nfix & mask(n) & eyedata(tr).stposrel(n)==eyedata(tr).stposrel(n+1) 
                            % multifix
                            k = k+1;
                            tmpfix.multifix  = k;
                            tmpfix.multi_n = multi_n;
                        else
                            if k>0
                                tmpfix.multifix  = k+1;
                                tmpfix.multi_n = multi_n;
                                multi_n = multi_n + 1;
                            else
                                tmpfix.multi_n = 0;
                                tmpfix.multifix  = k;
                            end
                            k = 0;
%                             if n<eyedata(tr).Nfix & mask(n) & eyedata(tr).stposrel(n)==eyedata(tr).stposrel(n-1) 
%                                 tmpfix.multifix  = k+1;
%                             else
%                                 tmpfix.multifix  = 1;
%                             end
                        end 
                        tmpfix.dur_added = tmpdur;
                        tmpfix.tonset  = eyedata(tr).fixs(n,1) - eyedata(tr).t_bgn_tr;%chequear
                        tmpfix.urevent = eyedata(tr).fixs(n,7);  
                        tmpfix.n       = n;
                        tmpfix.x       = eyedata(tr).fixs(n,4); 
                        tmpfix.y       = eyedata(tr).fixs(n,5);
                        tmpfix.istarget= eyedata(tr).fix2target(n); 
                        tmpfix.stposrel= eyedata(tr).stposrel(n);
                        tmpfix.stposabs  = eyedata(tr).stposabs(n);
                        tmpfix.stposrel= eyedata(tr).stposrel(n);
                        tmpfix.stposabs  = eyedata(tr).stposabs(n);
                        tmpfix.trial_targetseen     = any(eyedata(tr).fix2target);
                        tmpfix.trial_correct        = C;                      
                        tmpfix.trial_cat            = cat;
                        tmpfix.trial_resp           = resp;
                        tmpfix.trial_atr            = atr;
                        tmpfix.trial_type           = ExpTrials(tr).info.trialtype;
                        tmpfix.trial_cattype         = ExpTrials(tr).info.cattype;
                        tmpfix.trial_number         = tr;
                        tmpfix.trial_subject        = suj;
                        tmpfix.trial_bgn_latency    = eyedata(tr).t_bgn_tr;

                        if tmpfix.stposrel~=ExpTrials(tr).info.indtarget & tmpfix.stposabs~=0 ...
                                & (sum(ismember(ExpTrials(tr).info.indpos,tmpfix.stposabs))>0)

                            [cat_dist,atr_dist] = obj.fun_extract_dist_attributes(ExpTrials(tr), tmpfix.stposabs);
                                                    %esta funcion deberia arrojar cat y atr para distractor

                            tmpfix.isdistractor = 1;
                            tmpfix.distractor_cat          = cat_dist;
                            tmpfix.ditractor_atr           = atr_dist;
                        else
                            tmpfix.isdistractor = 0;
                            tmpfix.distractor_cat          = NaN;
                            tmpfix.ditractor_atr           = NaN;
                        end

                        fixs = [fixs, tmpfix];
                        
                    end
                    % Struct trials (all subjects)
                    if any(eyedata(tr).fix2target)
                        tmptrial.nfix2target    = find(eyedata(tr).fix2target,1,'first');
                        %analizar criterio (n al target tal que fixdur > a 
                        %el promedio de fixdur para distractores) )
                    else
                        tmptrial.nfix2target    = NaN;
                    end
                    tmptrial.Nfix           = eyedata(tr).Nfix;
                    tmptrial.targetseen     = any(eyedata(tr).fix2target);
                    tmptrial.correct        = C;
                    tmptrial.cat            = cat;
                    tmptrial.resp           = resp;
                    tmptrial.atr            = atr;
                    tmptrial.type           = ExpTrials(tr).info.trialtype;
                    tmptrial.number         = tr;
                    tmptrial.subject        = suj;

                    trials                  = [trials, tmptrial];

                    %CREO vector con las refijaciones para definir tmpfix.refix

                    ind         =  find([fixs.trial_number]==tr ...
                        & arrayfun( @(x)strcmp(x,suj),{fixs.trial_subject}));
                    tmp_pos_rel = [fixs(ind).stposrel];
                    tmp_refix   = zeros(size(tmp_pos_rel));

                    for i = 1:15
                        ind_refix = find(tmp_pos_rel==i);
                        if length(ind_refix)>1
                            tmp_refix(ind_refix) =  1:length(ind_refix);
                        end
                    end
                    ref              = num2cell([tmp_refix]);
                   [fixs(ind).refix] = deal(ref{:});    


                end

            end 
        end

    end
    
    methods(Static)
        %--------------------------------------------------------------------------------------------------------------------------------------------
        function [ha, pos]                              = tight_subplot(Nh, Nw, gap, marg_h, marg_w)  % mdf 14/04/21 to change space btw subplots
        % tight_subplot creates "subplot" axes with adjustable gaps and margins
        %
        % [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
        %
        %   in:  Nh      number of axes in hight (vertical direction)
        %        Nw      number of axes in width (horizontaldirection)
        %        gap     gaps between the axes in normalized units (0...1)
        %                   or [gap_h gap_w] for different gaps in height and width 
        %        marg_h  margins in height in normalized units (0...1)
        %                   or [lower upper] for different lower and upper margins 
        %        marg_w  margins in width in normalized units (0...1)
        %                   or [left right] for different left and right margins 
        %
        %  out:  ha     array of handles of the axes objects
        %                   starting from upper left corner, going row-wise as in
        %                   subplot
        %        pos    positions of the axes objects
        %
        %  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
        %           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
        %           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
        % Pekka Kumpulainen 21.5.2012   @tut.fi
        % Tampere University of Technology / Automation Science and Engineering
        if nargin<3; gap = .02; end
        if nargin<4 || isempty(marg_h); marg_h = .05; end
        if nargin<5; marg_w = .05; end
        if numel(gap)==1
            gap = [gap gap];
        end
        if numel(marg_w)==1
            marg_w = [marg_w marg_w];
        end
        if numel(marg_h)==1
            marg_h = [marg_h marg_h];
        end
        axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
        axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
        py = 1-marg_h(2)-axh; 
        % ha = zeros(Nh*Nw,1);
        ii = 0;
        for ih = 1:Nh
            px = marg_w(1);

            for ix = 1:Nw
                ii = ii+1;
                ha(ii) = axes('Units','normalized', ...
                    'Position',[px py axw axh], ...
                    'XTickLabel','', ...
                    'YTickLabel','');
                px = px+axw+gap(2);
            end
            py = py-axh-gap(1);
        end
        if nargout > 1
            pos = get(ha,'Position');
        end
        ha = ha(:);
        end
        function names                                  = loadSubjects(folder, expe, pre, type, keyw)     
            %Names of files to load (pre and expe are only used when data were recorded separately)
            % Load filenames (*.type) from folder.
            % if expe: loads suj.bdf/set 
            % if pre: loads suj_pre1/2.bdf/set
            % if keyw  = 'sth' loads files with 'sth' within file name
            names = dir(folder); 
            names = {names.name};

            if (nargin==5)
                %names = names(1,1:6);
                names = names(cellfun(@(x) ~isempty(x),regexp(names, keyw)));
            else
                names = dir([folder  '*.' type]);
                names = {names.name};
            end

            if pre == 1 && expe == 1
                return
            elseif pre == 1 
                names = names(cellfun(@(x) ~isempty(x),regexp(names, '_pre')));
            elseif expe == 1 
                names = names(cellfun(@(x) isempty(x),regexp(names, '_pre')));
            end  
        end
        function [vec_succesiv,ind_suc,mask]            = fun_detect_succesive_fixs(eyetrial)
            %    This function detect succesive fixations within the same trial. 
            %    
            %    input: 
            %          eyetrial: struct from eyedata(tr) 
            %
            %    output:
            %          vec_succesiv: vector         
            %          ind_suc:      indices  of succesive fixations
            %          mask:         boolean vector with succesive
            %                        fixations
            %    children:
            %            eyedata(tr)
            %Damian Care 27/6/19: 
            relpos=eyetrial.stposrel;
            %relpos= [1  0 2 2 2 7 7 0 0 2 2 9 0]';
            pri=[relpos~=0 & [diff(relpos);relpos(end)]==0] ;
            ult=[relpos~=0 & [relpos(1);diff(relpos)]==0];
            % mdf 3/6/21 vector with zero when fixs is not multifix and the position
            % of the fix when multifix>0
            vec_succesiv = relpos.*((pri+ult)>0);
            % mdf 3/6/21 mask with 1 for fixs of multifix, 0 otherwise
            mask = (pri+ult)>0;
            % mdf 3/6/21 indixes of succesive fixs in trial
            ind_suc = find((pri+ult)>0);
        end
        function [non_succesive_mask]                   = select_non_succesive_fixs(eyetrial,threshold)
            %    This function detect the non-succesive fixations within the same trial
            %    and select fixs greater than threshold
            %    
            %    input: 
            %          eyetrial: struct from eyedata(tr) 
            %          threshold: time 
            %
            %    output:
            %          non_succesive_mask: ones for fixs that are not
            %          multifix
            %    children:
            %            eyedata(tr)
            % Maria da Fonseca 7/6/21
            relpos=eyetrial.stposrel;
            %relpos= [1  0 2 2 2 7 7 0 0 2 2 9 0]';
            pri=[relpos~=0 & [diff(relpos);relpos(end)]==0] ;
            ult=[relpos~=0 & [relpos(1);diff(relpos)]==0];
            % mdf 3/6/21 vector with zero when fixs is not multifix and the position
            % of the fix when multifix>0
            vec_succesiv = relpos.*((pri+ult)>0);
            % mdf 3/6/21 mask with 1 for fixs of multifix, 0 otherwise
            mask = (pri+ult)>0;
            
            % array with fix duration
            dur = []; 
            Ntr = length(relpos);
            for n = 1:Ntr
                dur= [dur eyetrial.fixs(n,2) - eyetrial.fixs(n,1)]; 
            end
            mask_dur = dur>threshold;
            
            non_succesive_mask = (~mask)'.*mask_dur;
        end
        function [fix2keep,ind_suc,mask]                = select_multifix(eyetrial,threshold)
            %    This function detect succesive fixations within the same trial
            %    and select the first fix greater than threshold
            %    
            %    input: 
            %          eyetrial: struct from eyedata(tr) 
            %          threshold: time 
            %
            %    output:
            %          fix2keep: position of the fixs to keep
            %          ind_suc: index of fixs to keep 
            %          mask: ones for fix to keep, zero otherwise 
            %
            %    children:
            %            eyedata(tr)
            % Maria da Fonseca 7/6/21
            relpos=eyetrial.stposrel;
            %relpos= [1  0 2 2 2 7 7 0 0 2 2 9 0]';
            pri=[relpos~=0 & [diff(relpos);relpos(end)]==0] ;
            ult=[relpos~=0 & [relpos(1);diff(relpos)]==0];
            % mdf 3/6/21 vector with zero when fixs is not multifix and the position
            % of the fix when multifix>0
            vec_succesiv = relpos.*((pri+ult)>0);
            % mdf 3/6/21 mask with 1 for fixs of multifix, 0 otherwise
            mask = (pri+ult)>0;
            % array with fix duration
            dur = []; 
            Ntr = length(relpos);
            for n = 1:Ntr
                dur= [dur eyetrial.fixs(n,2) - eyetrial.fixs(n,1)]; 
            end
            mask_dur = dur>threshold;
            
            x = vec_succesiv'.*mask_dur;
            fix2keep = x.*((x~=0) & diff([0 x]));
            mask = fix2keep>0;
            ind_suc = find(fix2keep>0);
        end
        function [categoria,atributo,tipo_respuesta]    = fun_extract_attributes_target(SingleTrial)
            %    This function retrieve stimulus data coming from a target. 
            %    
            %    input: 
            %          eyetrial:      SingleTrial struct from Exp
            %
            %    output:
            %          categoria:     Category of the target
            %          atributo:      Target atribute
            %          tipo_respuesta:P for present stimulus and A for
            %          absent.
            %
            %    
            %Damian Care 27/6/19: 
            %disp(SingleTrial.info.filenames_target)
            try
            s = regexp(SingleTrial.info.filenames_target, ['../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+-\w+.jpg|'...
                '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+.jpg|'...
                '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+.\w+.jpg|'...
                '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+-\w+-\w+.jpg|'...
                '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+-\w+-\w+-.jpg'], 'names');
            
                categoria   = s.categoria;
                atributo    = s.atributo; 
            catch
                disp(s)
                keyboard
            end

            if any(ismember(SingleTrial.info.filenames_grid,SingleTrial.info.filenames_target))
                tipo_respuesta = 'P';
            else
                tipo_respuesta = 'A';
%                 categoria = [];%agregado para lidiar con error J03 y J06
%                 atributo  = [];%" "            "             "        "
            end
        end
        function [cat_dist,atr_dist]                    = fun_extract_dist_attributes(SingleTrial, abs_position)
            %    This function retrieve stimulus data coming from the ditractor
            %    at the given position.
            %    SingleTrial data should be a matlab struct.
            %
            %    input: 
            %          SingleTrial :      SingleTrial struct from ExpTrials
            %          rel_position:      Number of the distractor on the grid  
            %                             of stimuli.        
            %
            %    output:
            %          cat_dist:     Category of the distractor
            %          atr_dist:     Distractor atribute
            %      
            %
            %    
            %Damian Care 30/1/21: 
            %disp(SingleTrial.info.filenames_target)
            %disp(SingleTrial.info.filenames_target)



            dist_ind = find(ismember(SingleTrial.info.indpos,abs_position));

          try
            s = regexp(SingleTrial.info.filenames_grid{dist_ind}, ['../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+-\w+.jpg|'...
                '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+.jpg|'...
                '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+.\w+.jpg|'...
                 '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+-\w+-\w+-.jpg|'...
                '../experiment/(?<categoria>\w+)/(?<atributo>\w+)_shined_luminance2/\w+-\w+-\w+.jpg'], 'names');
            %keyboard
                cat_dist   = s.categoria;
                atr_dist    = s.atributo; 
          catch
              fprintf(SingleTrial.info.filenames_grid{dist_ind})
          end
        end
        function eyedata                                = fun_parse_eyedata_EEG_beta(EEG,ExpTrials)
            % Version of the fun_parse_eyedata function using eyetracking data from
            % EEG struct of EEGLAB. It creates a struct with datapaths for
            % each trial, position of stimuli number of total fixations and
            % fixations to target.
            %
            %input:
            %      EEG: EEGLAB struct.
            %      ExpTrials: struct from behavioral data from experiment.
            %
            %output: 
            %       eyedata: struct with eyetracker data arrange acording
            %       trials. Infomation about number of fixations and
            %       stimuli positions are also reported.
            % 
            %
            % 2021-21-01 Damian Care
            % Eye Movements
            % Then I can check for sustained fixations during Target presentation, 
            % fixation period, and retain period, or eye movements in the eyemap
            % periods

            padding_pre     = 0;%EEG.xmin;
            padding_post    = 0;

            event = EEG.event;

            indini = find(ismember({event.type},{'13'}));
            %indend = find();
        %     trialdur = [[event([indend]).latency]'-[event([indini]).latency]']/EEG.srate;    

            t_bgn = [event( indini).latency]'/EEG.srate;
            t_end = [event( indini).latency]'/EEG.srate + EEG.xmax;


            eyedata = [];


            indfix = ismember({event.type},{'fixation'});
            %find fixations inside the trial
        %   indsac = ismember({event.type},{'saccade'});

            fixEEG = [[event(indfix).latency]'/EEG.srate,...                            %onset
                     ([event(indfix).latency]'+ [event(indfix).duration]')/EEG.srate,...%offset
                     [event(indfix).duration]'/EEG.srate,...                            %duration
                     [event(indfix).fix_avgpos_x]',[event(indfix).fix_avgpos_y]',...    %position x and y
                     [event(indfix).fix_avgpupilsize ]',...                             %pupil
                     [event(indfix).urevent ]'];                                        %urevent
                 


            if length(ExpTrials) == sum(ismember({EEG.event.type},{'13'}))
                for tr = 1:length(ExpTrials)
                    indtr = [];
                    info = ExpTrials(tr).info;         

                    % fix: t_start t_end dur x_medio y_medio pupila_medio
                    indfix_tr = [event(indfix).epoch]==tr & strcmp({event(indfix).type},'fixation')... belongs to trial tr and is a fixation
                        & t_bgn(tr)*EEG.srate + padding_pre  <[event(indfix).latency] &... %fall within trial 
                        t_end(tr)*EEG.srate>[event(indfix).latency]; 
                    Nfix      = sum(indfix_tr);
                    if (Nfix>0)
                        fixations   = nan(Nfix, size(fixEEG,2));
                        for i = 1:sum(indfix_tr)
                            indtr = find(indfix_tr);
                            fixations(i,:)  = fixEEG(indtr(i),:);
                            %keyboard
                        end
            % Complete for saccades
                        
                        xpos = fixations(:,4);
                        ypos = fixations(:,5);
                        stposrel = zeros(Nfix,1);
                        stposabs = zeros(Nfix,1);
                        for j = 1:Nfix
                            for i = 1:length(info.dstRects)
                                if IsInRect(xpos(j),ypos(j),info.dstRects(:,i))  
                                    stposrel(j) = i;%info.indpos(i); % mdf 17/6/21 index position
                                    stposabs(j) = info.indpos(i); % mdf 17/6/21 position number, 0 for non-stims places
                                    break;
                                end
                            end
                        end

                        if (info.indtarget>0)
                            fix2target = zeros(Nfix,1);
                            for j = 1:Nfix
                                fix2target(j) = IsInRect(xpos(j),ypos(j),info.recttarget);
                            end
                        else
                            fix2target = zeros(Nfix,1);
                        end

                        eyedata(tr).Nfix      = Nfix;
                        eyedata(tr).fixs      = fixations;
                        %eyedata(tr).Nsac      = Nsac;
                        %eyedata(tr).sacs      = saccades;
                        eyedata(tr).stposrel  = stposrel;
                        eyedata(tr).stposabs  = stposabs;
         


                        eyedata(tr).fix2target= fix2target;
                        eyedata(tr).t_bgn_tr = t_bgn(tr);
                        %eyedata(tr).samples   = all.samples(all.samples(:,1)>t_bgn(tr) & all.samples(:,1)<t_end(tr),:);
                        %eyedata(tr).samples(:,1) = eyedata(tr).samples(:,1)-eyedata(tr).samples(1,1);

                        % AGREGAR DECIMATION Y DECLARAR LA NUEVA sampling rate
                        % eyedata(tr).srate       = ????
                    else
                        %indsac = (tmpsac(:,2)>t_bgn(tr) & tmpsac(:,1)<t_end(tr));
                        %Nsac = sum(indsac);

                        eyedata(tr).Nfix      = Nfix;
                        eyedata(tr).fixs      = [];
                        %eyedata(tr).Nsac      = Nsac;
                        %eyedata(tr).sacs      = [];
                        eyedata(tr).stposrel  = [];
                        eyedata(tr).stposabs  = [];
                        eyedata(tr).fix2target= [];
                        %eyedata(tr).samples   = [];
                        eyedata(tr).urevent    = [event(find(indfix)).urevent];%para agregar luego info al EEG
                    end
                    eyedata(tr).padding_pre     = padding_pre;
                    eyedata(tr).padding_post    = padding_post;
            end
            else
                fprintf('number of trials does not match btw mat and EEG files')
        end
        end
        function ERPstruct                              = fun_makeStructToERP(ERPstruct,condVec,iSuj,folder,subjName)
            % This function appends ERPs from a subject to an array of ERPs.
            % condition epochs should be previously created using indexes
            % from indexCondEpochs5(EEG,fixs,  durMinFix) function
            % 
            % 
            %
            %input:
            %      EEGstruct: [] (empty list or list of existing ERPs to be append to)
            %      condVec: cell array with condition names of each ERP (e.g.: 
            %      condVec = {'_NtEXO.set' '_NtEXF.set' '_NtEXOsinI.set' '_NtEXFsinI.set'})
            %      iSuj: number identifying subject to be consider for ERP
            %      calculation.
            %      folder: path where condition epochs are stored.
            %      subjName: name of the subjet to be stored in the output
            %      struct.
            %
            %output: 
            %       ERPstruct: ERPstruct with ERPs from iSuj appended.
            %       
            %       
            % 
            %
            % 2021-21-01 Damian Care
                filename    = [folder subjName '_fixEpoch.set'];    

                EEG     = pop_loadset(filename);
                EEG     = eeg_checkset( EEG ); 
                ERPstruct(iSuj).times = EEG.times;
                ERPstruct(iSuj).chanlocs = EEG.chanlocs(1:64);
                ERPstruct(iSuj).Nepoch = size(EEG.data,3);
                ERPstruct(iSuj).fixepoch = mean(EEG.data(1:64,:,:),3);
                ERPstruct(iSuj).Suj      = subjName;

                %Distractores Caras
                    
                filename    = [folder subjName condVec{1}];
                if exist(filename)
                EEG1     = pop_loadset(filename);
                EEG1     = eeg_checkset( EEG1 );
                ERPstruct(iSuj).N1 = size(EEG1.data,3);
                ERPstruct(iSuj).(condVec{1}(2:end-4)) = mean(EEG1.data(1:64,:,:),3);

                end
                %Distractores Objetos
                filename    = [folder subjName condVec{2}];
                if exist(filename)
                EEG2     = pop_loadset(filename);
                EEG2     = eeg_checkset( EEG2 );
                ERPstruct(iSuj).N2 = size(EEG2.data,3);
                ERPstruct(iSuj).(condVec{2}(2:end-4)) = mean(EEG2.data(1:64,:,:),3);

                end
                %Target Objetos
                filename    = [folder subjName  condVec{3}];
                if exist(filename)                
                EEG3     = pop_loadset(filename);
                EEG3     = eeg_checkset( EEG3 );
                ERPstruct(iSuj).N3 = size(EEG3.data,3);
                ERPstruct(iSuj).(condVec{3}(2:end-4)) = mean(EEG3.data(1:64,:,:),3);

                end
                %Target Caras
                filename    = [folder subjName  condVec{4}];
                if exist(filename)
                EEG4     = pop_loadset(filename);
                EEG4     = eeg_checkset( EEG4 );
                ERPstruct(iSuj).N4 = size(EEG4.data,3);
                ERPstruct(iSuj).(condVec{4}(2:end-4)) = mean(EEG4.data(1:64,:,:),3);

                end
                % suj x cond x canales x tiempo
               
           
        end
        function ERPstruct                              = fun_makeStructToERP2021(ERPstruct,cond1,cond2,iSuj,folder,subjName)
            % This function appends ERPs from a subject to an array of ERPs.
            % condition epochs should be previously created using indexes
            % from indexCondEpochs5(EEG,fixs,  durMinFix) function
            % (This is a 2021 version of the previous one that creates only one 
            % struct for condition to analyse).
            % 
            % 
            %
            %input:
            %      EEGstruct: [] (empty list or list of existing ERPs to be append to)
            %      cond: string with condition names of each ERP (e.g.: 
            %      cond = 'NTO.set' it should be the name asociated to EEG
            %      epoch data.
            %      iSuj: number identifying subject to be consider for ERP
            %      calculation.
            %      folder: path where condition epochs are stored.
            %      subjName: name of the subjet to be stored in the output
            %      struct.
            %
            %output: 
            %       ERPstruct: ERPstruct with ERPs from iSuj appended.
            %       
            %       
            % 
            %
            % 2021-21-01 Damian Care
                
                filename    = [folder subjName '_fixEpoch.set'];    

                EEG     = pop_loadset(filename);
                EEG     = eeg_checkset( EEG ); 
                ERPstruct(iSuj).times = EEG.times;
                ERPstruct(iSuj).chanlocs = EEG.chanlocs(1:64);
                ERPstruct(iSuj).n_fixsepoch = size(EEG.data,3);
                ERPstruct(iSuj).fixepoch = mean(EEG.data(1:64,:,:),3);
                ERPstruct(iSuj).Suj      = subjName;

                %Distractores Caras
                    
                filename    = [folder subjName cond1];

                if exist(filename)
                EEG1     = pop_loadset(filename);
                EEG1     = eeg_checkset( EEG1 );

                ERPstruct(iSuj).(['n_' (cond1(2:end-4))]) = size(EEG1.data,3);
                ERPstruct(iSuj).(cond1(2:end-4)) = mean(EEG1.data(1:64,:,:),3);
               
                end
                %Distractores Objetos
                filename    = [folder subjName cond2];
                if exist(filename)
                EEG2     = pop_loadset(filename);
                EEG2     = eeg_checkset( EEG2 );
                ERPstruct(iSuj).(['n_' (cond2(2:end-4))]) = size(EEG2.data,3);
                ERPstruct(iSuj).(cond2(2:end-4)) = mean(EEG2.data(1:64,:,:),3);
                
                end           
        end
        function ERPstruct                              = fun_makeStructToERP2021_extended(ERPstruct,c1,c2,c3,c4,c5,c6,iSuj,folder,subjName)
            % This function appends ERPs from a subject to an array of ERPs.
            % condition epochs should be previously created using indexes
            % from indexCondEpochs5(EEG,fixs,  durMinFix) function
            % (This is a 2021 version of the previous one that creates only one 
            % struct for condition to analyse).
            % 
            % 
            %
            %input:
            %      EEGstruct: [] (empty list or list of existing ERPs to be append to)
            %      cond: string with condition names of each ERP (e.g.: 
            %      cond = 'NTO.set' it should be the name asociated to EEG
            %      epoch data.
            %      iSuj: number identifying subject to be consider for ERP
            %      calculation.
            %      folder: path where condition epochs are stored.
            %      subjName: name of the subjet to be stored in the output
            %      struct.
            %
            %output: 
            %       ERPstruct: ERPstruct with ERPs from iSuj appended.
            %       
            %       
            % 
            %
            % 2021-21-01 Damian Care
                
                filename    = [folder subjName '_fixEpoch.set'];    

                EEG     = pop_loadset(filename);
                EEG     = eeg_checkset( EEG ); 
                ERPstruct(iSuj).times = EEG.times;
                ERPstruct(iSuj).chanlocs = EEG.chanlocs(1:64);
                ERPstruct(iSuj).n_fixsepoch = size(EEG.data,3);
                ERPstruct(iSuj).fixepoch = mean(EEG.data(1:64,:,:),3);
                ERPstruct(iSuj).Suj      = subjName;

                %c1
                filename    = [folder subjName c1];
                if exist(filename)
                    EEG1     = pop_loadset(filename);
                    EEG1     = eeg_checkset( EEG1 );
                    ERPstruct(iSuj).(['n_' (c1(2:end-4))]) = size(EEG1.data,3);
                    ERPstruct(iSuj).(c1(2:end-4)) = mean(EEG1.data(1:64,:,:),3);
                end
                
                %c2
                filename    = [folder subjName c2];
                if exist(filename)
                    EEG2     = pop_loadset(filename);
                    EEG2     = eeg_checkset( EEG2 );
                    ERPstruct(iSuj).(['n_' (c2(2:end-4))]) = size(EEG2.data,3);
                    ERPstruct(iSuj).(c2(2:end-4)) = mean(EEG2.data(1:64,:,:),3);
                end  
                
                %c3
                filename    = [folder subjName c3];
                if exist(filename)
                    EEG3     = pop_loadset(filename);
                    EEG3     = eeg_checkset( EEG3 );
                    ERPstruct(iSuj).(['n_' (c3(2:end-4))]) = size(EEG3.data,3);
                    ERPstruct(iSuj).(c3(2:end-4)) = mean(EEG3.data(1:64,:,:),3);
                end
                
                %c4
                filename    = [folder subjName c4];
                if exist(filename)
                    EEG4     = pop_loadset(filename);
                    EEG4     = eeg_checkset( EEG4 );
                    ERPstruct(iSuj).(['n_' (c4(2:end-4))]) = size(EEG4.data,3);
                    ERPstruct(iSuj).(c4(2:end-4)) = mean(EEG4.data(1:64,:,:),3);
                end
                
                %c5
                filename    = [folder subjName c5];
                if exist(filename)
                    EEG5     = pop_loadset(filename);
                    EEG5     = eeg_checkset( EEG5 );
                    ERPstruct(iSuj).(['n_' (c5(2:end-4))]) = size(EEG5.data,3);
                    ERPstruct(iSuj).(c5(2:end-4)) = mean(EEG5.data(1:64,:,:),3);
                end
                
                %c6
                filename    = [folder subjName c6];
                if exist(filename)
                    EEG6     = pop_loadset(filename);
                    EEG6     = eeg_checkset( EEG6 );
                    ERPstruct(iSuj).(['n_' (c6(2:end-4))]) = size(EEG6.data,3);
                    ERPstruct(iSuj).(c6(2:end-4)) = mean(EEG6.data(1:64,:,:),3);
                end
        end
        function [index,EEG]                            = indexCondEpochs(EEG,freqFilteredEEG, fixsSuj,minfixdur,limits,args) 
            % Function that produce lists of indexes to be consider
            % according to the disired conditions defined----
            % 
            % 
            % 
            %
            %input:
            %      EEG: EEGlab with epochs from all fixations of one
            %      freqFilteredEEG: 1 if the EEG is filtered in freq
            %      subject.
            %      fixsSuj:   Struct containing all matadata for each fixation
            %      minfixdur: minimun duration cut for fixations to be consider
            %      subjName:  name of the subjet to be stored in the output
            %      struct.
            %      exp:       VS,  EX. If empty both are consider.(string)
            %      trial:     'O' = object ,'F' = faces upward , 'I' = inverted faces,(string)
            %      stim:      df = dist face , do = dist obj,tf = target face, to= tar obj(string).
            %      correct:   0 = incorret answer,1 = correct answer. If
            %      absent both are consider. (int)
            %      present:   'A' = target absent,'P' = target present. If
            %      absent both are consider. (string)
            %
            %output: 
            %       index:    indexes of epoch for selected conditions.
            %       EEG:      eeglab struct as in input.
            %       
            %children:
            %         fixationsAnalysisOne(): EEG and fixsSuj are outputs
            %         of this function.
            % 
            %future: so far I have a different function for different
            %conditions I plan to add new variables to use only one
            %function.
            % 2021-21-01 Damian Care
            event = EEG.event;
            fixs = fixsSuj;

            EEG = pop_epoch( EEG, {  'fixation'  }, limits, 'newname', 'fixs and conditions epoch', 'epochinfo', 'yes');%epoqueo por fijaciones para el sujeto suj
            if freqFilteredEEG
                EEG = pop_rmbase( EEG, [-199.2 0]);
            else
                EEG = pop_rmbase( EEG, [-199.2 0]);%baseline to 0 secs 27/02/2020 
            end 
            
            EEG.trials
            % 
            %keyboard
            %urevents = arrayfun(@(y) y.eventurevent{cellfun(@(x) x==0, y.eventlatency)},EEG.epoch); %busco eventos de la fijacion principal
            %urevents  = arrayfun(@(y) y.eventurevent{cellfun(@(x) strcmp(x,'fixation'), y.eventtype) & cellfun(@(x) x==0, y.eventlatency)},EEG.epoch);
            %urevents1 = arrayfun(@(y) y.eventurevent{cellfun(@(x) strcmp(x,'fixation'), y.eventtype)},EEG.epoch);
            %urevents2 = arrayfun(@(y) y.eventurevent{cellfun(@(x) x==0, y.eventlatency)},EEG.epoch);
            
            urevent_fix = [];
            urevent_lat = [];
            for i=1:EEG.trials
                urevent_fix = [urevent_fix EEG.epoch(i).eventurevent{cellfun(@(x) strcmp(x,'fixation'), EEG.epoch(i).eventtype)}];
                urevent_lat = [urevent_lat EEG.epoch(i).eventurevent{cellfun(@(x) x==0, EEG.epoch(i).eventlatency)}];
            end
            urevents = intersect(urevent_fix,urevent_lat);
            
%             fprintf('num of fixations %i, num of events with lat 0 %i, num of intersection %i \n', numel(urevent_fix),...
%                 numel(urevent_lat),numel(urevents));
            
            infix = ismember([fixs.urevent],urevents); % indices de las fijaciones que corresponden al fixs con urevents del epoqueado por fixs
%             fprintf('num of infix %i \n', numel(infix));
            
            %index to retain
            index     = [];
                        
            
            stimsdist = {'NTF' 'NTO'};
            stimstarg = {'TF' 'TO'};
            
            exps     = {'VS', 'EX'};  
            trials = {'F','O','I'};
            corrects = [1 0];      
            presents = {'P','A'};
            expbypass = 0;
            trialbypass = 0;
            correctbypass = 0;
            presentbypass = 0;
            trialcatbypass = 0; % mdf 3/5/21 added for trials in which participats were looking for F and I faces
            tonsetbypass = 0;  % mdf 3/5/21 added for limit fixations with t_onset > args.tonsetlim
            multibypass = 0; % mdf 10/5/21 added to avoid fix with multifix>1
            pretargetbypass = 0; % mdf 18/6/21 to filter fixs previous to target
            
             if ~isfield(args,'fixorder')
                 fixbypass = 1;
                 fixorder = 0;
             else
                 fixbypass = 0;
                 fixorder = args.fixorder;
             end       
            
             if isfield(args,'stim')
                stim = args.stim;
                distmask = ismember(stimsdist, stim);
                targmask = ismember(stimstarg, stim);
                if sum(distmask)==1
                    isdist   = 1;
                    istarg = 0;
                    targcat = [];
                    targcatbypass = 1;
                    distcatbypass = 0;
                    if distmask(1) == 1
                      distcat = 'faces';
                    elseif distmask(2) == 1
                      distcat = 'objects';
                    else
                      distcat = [];
                    end
                elseif sum(targmask) ==1
                    isdist   = 0;
                    istarg= 1;
                    distcat = [];
                    distcatbypass = 1;
                    targcatbypass = 0;
                    if targmask(1) == 1
                      targcat = 'faces';
                    elseif targmask(2) == 1
                      targcat = 'objects';
                    else
                      targcat = [];
                    end
                elseif sum(distmask) == 2
                    isdist   = 1;
                    istarg = 0;
                    distcat = [];
                    targcat = [];
                    targcatbypass = 1;
                    distcatbypass = 1;
                end
                stimbypass = 0;

             else
                 stimbypass = 1;
                 targcatbypass = 1;
                 distcatbypass = 1;
                 stim = [];
                 isdist = 0;
                 istarg = 0;
                 distcat = [];
                 targcat = [];
             end
             
        
             if ~isfield(args,'exp')
                expbypass = 1;
                exp = [];
             else
                 exp = args.exp;
                 expbypass = 0;
             end
             
             if ~isfield(args,'trial_cat') % mdf 3/5/21 now is only considering faces, for objects we can use args.trial='O'
                trialcatbypass = 1;
                trial_cat = 0;
             else
                trial_cat = args.trial_cat;
             end

             if ~isfield(args,'trial')
                trialbypass = 1;
                trial = 0;
             else
                trial = args.trial;
             end

             
             if ~isfield(args,'correct')
                correctbypass = 1;
                correct = 0;
             else
                 correct = args.correct;
                 correctbypass = 0;
             end
             
             if ~isfield(args,'present')
                presentbypass = 1;
                present = 0;
             else
                 present = args.present;
                 presentbypass = 0;
             end    
             
%              if ~isfield(args,'rank')
%                  rankbypass = 1;
%                  rank = 0;
%              else
%                  rank = args.rank;
%                  rankbypass = 0;
%              end    
             
             if ~isfield(args,'tonsetlim')
                 tonsetbypass = 1;
                 tonset = 0;
             else
                 tonset = args.tonsetlim;
             end
             
             if ~isfield(args,'multifix')
                 multibypass = 1;
                 multi = 0;
             else
                 multi = args.multifix;
             end
             
             if ~isfield(args,'pretarget')
                pretargetbypass = 1;
                pretarget = 0;
             else
                pretargetbypass = 0;
                pretarget = args.pretarget;
             end

            tmax = minfixdur ;
            indfixnum =  find(infix);
            eventos = [];
            for i = 1:sum(infix)
                    in = indfixnum(i);
                    urin   = find([urevents]== fixs(in).urevent);
                    if ((fixs(in).refix==0 || fixs(in).refix==1) && ...
                        (fixs(in).dur_added>=tmax) &&...
                        (fixs(in).rank>1) &&...
                        (fixs(in).isdistractor == isdist || stimbypass )&&...
                        (strcmp(fixs(in).distractor_cat, distcat) || stimbypass ||distcatbypass) &&...
                        (fixs(in).istarget == istarg || stimbypass) && ...
                        (strcmp(fixs(in).trial_type, exp) || expbypass ) &&... 
                        (strcmp(fixs(in).trial_cattype, trial)|| trialbypass || strcmp(fixs(in).trial_cat, trial)) &&...
                        (strcmp(fixs(in).trial_cat, targcat) || stimbypass || targcatbypass) &&...
                        (fixs(in).trial_correct == correct || correctbypass) &&...
                        (fixs(in).trial_resp ==present || presentbypass) &&...
                        (strcmp(fixs(in).trial_cat, trial_cat) || trialcatbypass) &&...
                        (fixs(in).tonset < tonset || tonsetbypass) &&...
                        (fixs(in).multifix == multi || multibypass) &&...
                        (fixs(in).rank == fixorder || fixbypass) &&...
                        (fixs(in).pretarget == pretarget || pretargetbypass))
                        %keyboard
                        %indexs to retain
                        index = [index urin];
                        eventos = [eventos fixs(in).urevent];
                        
                    end
            end
        end
        function [index,EEG]                            = index_from_urevents(EEG, fix_urevents,minfixdur,limits) 
            EEG = pop_epoch( EEG, {  'fixation'  }, limits, 'newname', 'fixs and conditions epoch', 'epochinfo', 'yes');%epoqueo por fijaciones para el sujeto suj
            EEG = pop_rmbase( EEG, [-199.2188 0]);%baseline to 0 secs 27/02/2020
            
            urevent_fix = [];
            urevent_lat = [];
            for i=1:EEG.trials
                urevent_fix = [urevent_fix EEG.epoch(i).eventurevent{cellfun(@(x) strcmp(x,'fixation'), EEG.epoch(i).eventtype)}];
                urevent_lat = [urevent_lat EEG.epoch(i).eventurevent{cellfun(@(x) x==0, EEG.epoch(i).eventlatency)}];
            end
            urevents = intersect(urevent_fix,urevent_lat);
            
            infix = ismember(urevents,fix_urevents); % indices de las fijaciones que corresponden al fixs con urevents del epoqueado por fixs  

            index = find(infix==1);
        end        
        function [index,EEG]                            = prepareEEG2LMM(EEG, fixsSuj,minfixdur,args) 
            % Function that produce lists of indexes to be
            % considerXXXXXXWRONG
            % according to the disired conditions defined----
            % 
            % 
            % 
            %
            %input:
            %      EEG: EEGlab with epochs from all fixations of one
            %      subject.
            %      fixsSuj:   Struct containing all matadata for each fixation
            %      minfixdur: minimun duration cut for fixations to be consider
            %      subjName:  name of the subjet to be stored in the output
            %      struct.
            %      exp:       VS,  EX. If empty both are consider.(string)
            %      trial:     'O' = object ,'F' = faces upward , 'I' = inverted faces,(string)
            %      stim:      df = dist face , do = dist obj,tf = target face, to= tar obj(string).
            %      correct:   0 = incorret answer,1 = correct answer. If
            %      absent both are consider. (int)
            %      present:   'A' = target absent,'P' = target present. If
            %      absent both are consider. (string)
            %
            %output: 
            %       index:    indexes of epoch for selected conditions.
            %       EEG:      eeglab struct as in input.
            %       
            %children:
            %         fixationsAnalysisOne(): EEG and fixsSuj are outputs
            %         of this function.
            % 
            %future: so far I have a different function for different
            %conditions I plan to add new variables to use only one
            %function.
            % 2021-21-01 Damian Care
            event = EEG.event;
            fixs = fixsSuj;

% MDF 26/04/21 OJO QUE ACA HAY UN EPOCHEADO ENTRE -200 Y 800 MS...AHORA
% CAMBIAMOS A LA VENTANA -200,400 MS. 

            EEG = pop_epoch( EEG, {  'fixation'  }, [-0.2 0.8], 'newname', 'fixs and conditions epoch', 'epochinfo', 'yes');%epoqueo por fijaciones para el sujeto suj
            EEG = pop_rmbase( EEG, [-199.2188 0]);%baseline to 0 secs 27/02/2020 
            % 
            %keyboard
            urevents = arrayfun(@(y) y.eventurevent{cellfun(@(x) x==0, y.eventlatency)},EEG.epoch); %busco eventos de la fijacion principal
            infix = ismember([fixs.urevent],urevents); % indices de las fijaciones que corresponden al fixs con urevents del epoqueado por fixs



            %index to retain
            index     = [];
                        
            
            stimsdist = {'NTF' 'NTO'};
            stimstarg = {'TF' 'TO'};
            
            exps     = {'VS', 'EX'};  
            trials = {'F','O','I'};
            corrects = [1 0];      
            presents = {'P','A'};
            expbypass = 0;
            trialbypass = 0;
            correctbypass = 0;
            presentbypass = 0;
             if isfield(args,'stim')
                stim = args.stim;
                distmask = ismember(stimsdist, stim);
                targmask = ismember(stimstarg, stim);
                if sum(distmask)==1
                    isdist   = 1;
                    istarg = 0;
                    targcat = [];
                    targcatbypass = 1;
                    distcatbypass = 0;
                    if distmask(1) == 1
                      distcat = 'faces';
                    elseif distmask(2) == 1
                      distcat = 'objects';
                    else
                      distcat = [];
                    end
                elseif sum(targmask) ==1
                    isdist   = 0;
                    istarg= 1;
                    distcat = [];
                    distcatbypass = 1;
                    targcatbypass = 0;
                    if targmask(1) == 1
                      targcat = 'faces';
                    elseif targmask(2) == 1
                      targcat = 'objects';
                    else
                      targcat = [];
                    end
                end
                 stimbypass = 0;
             else
                 stimbypass = 1;
                 targcatbypass = 1;
                 distcatbypass = 1;
                 stim = [];
                 isdist = 0;
                 istarg = 0;
                 distcat = [];
                 targcat = [];
             end
             
        
             if ~isfield(args,'exp')
                expbypass = 1;
                exp = [];
             else
                 exp = args.exp;
                 expbypass = 0;
             end
             
             if ~isfield(args,'trial')
                trialbypass = 1;
                trial = 0;
             else
                trial = args.trial;
             end
             
             if ~isfield(args,'correct')
                correctbypass = 1;
                correct = 0;
             else
                 correct = args.correct;
                 correctbypass = 0;
             end
             
             if ~isfield(args,'present')
                presentbypass = 1;
                present = 0;
             else
                 present = args.present;
                 presentbypass = 0;
             end    
             
             if ~isfield(args,'rank')
                 rankbypass = 1;
                 rank = 0;
             else
                 rank = args.rank;
                 rankbypass = 0;
             end    
            
            tmax = minfixdur ;
            indfixnum =  find(infix);
            for i = 1:sum(infix)
                    in = indfixnum(i);
                    urin   = find([urevents]== fixs(in).urevent);
                    if ((fixs(in).refix==0 || fixs(in).refix==1) && ...
                        fixs(in).dur_added>=tmax &&...
                        (fixs(in).isdistractor == 1 || fixs(in).istarget == 1 )&&...%hardoded to include dist and targ
                        (strcmp(fixs(in).distractor_cat, distcat) || stimbypass ||distcatbypass) &&...
                        (strcmp(fixs(in).trial_type, exp) || expbypass ) &&... 
                        (strcmp(fixs(in).trial_cattype, trial)|| trialbypass) &&... 
                        (strcmp(fixs(in).trial_cat, targcat) || stimbypass || targcatbypass) &&...
                        (fixs(in).trial_correct == correct || correctbypass) &&...
                        (fixs(in).trial_resp ==present || presentbypass) &&...
                        (fixs(in).n == rank || rankbypass)) 
                        
                        EEG.epoch(urin).rank     = fixs(in).n;
                        EEG.epoch(urin).fixDur   = fixs(in).dur_added;
                        EEG.epoch(urin).trial    = fixs(in).trial_cattype;
                        EEG.epoch(urin).isDist   = fixs(in).isdistractor;
                        if fixs(in).isdistractor
                            EEG.epoch(urin).category  = fixs(in).distractor_cat;
                        else
                            EEG.epoch(urin).category  = fixs(in).trial_cat;
                        end
                        EEG.epoch(urin).taskType = fixs(in).trial_type;

                        
                        %keyboard
                        %indexs to retain
                        index = [index urin];
                    end
            %keyboard

            end
            EEG = pop_select(EEG,'trial',  index);
            EEG = pop_select(EEG,'time',[-0.2 EEG.xmax]);
            EEG = pop_select(EEG,'notrial', find(arrayfun(@(x) any(strcmp(x.eventtype,'bad_ET')),EEG.epoch)) );
        end
        function Results                                = ept_TFCE(DataFile1, DataFile2, ElecFile, varargin)
                %This is an imported version of the main TFCE function from  Mensen toolbox    
                %% EEG Permutation Test (egt) using Threshold-Free Cluster Enhancement 

                % Copyright(C) 2012  Armand Mensen (14.12.2010)

                % This program is free software: you can redistribute it and/or modify
                % it under the terms of the GNU General Public License as published by
                % the Free Software Foundation, either version 3 of the License, or
                % (at your option) any later version.
                % This program is distributed in the hope that it will be useful,
                % but WITHOUT ANY WARRANTY; without even the implied warranty of
                % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
                % GNU General Public License for more details.
                % You should have received a copy of the GNU General Public License
                % along with this program.  If not, see <http://www.gnu.org/licenses/>.

                % [Description]
                % This tests initially computes a T-Value for each channel and sample
                % These values are then enhanced or reduced depending on the size of
                % the T-value and the size of the cluster the channel belongs to at 
                % various thresholds...
                %     
                % TFCE is essentially the sum of the clusters to the power of 0.5 multiplied
                % by the height of the current threshold squared
                % 
                % [Input]
                % Participant data should be organised into two factors
                %     e.g. Dataset 1 is controls while Dataset 2 is patients
                % Electrode locations file created using eeglab is required to 
                % % calculate channel neighbours
                % 
                % Analysis Types
                % i = independent sample T-Test
                % d = dependent (paired) sample T-Test
                % o = one-sample T-test
                % c = Pearson's correlation (r)
                % 
                % This Summary file should be a "Participants x Channels x Samples" variable
                % - Participants should be in order of the group they belong to...
                % - Channels must be in the same order as the corresponding electrodes file
                % - Samples should be in chronological order
                % - Values in the variable correspond to the mV averaged ERP data
                % 
                % [Output]
                % - Info Structure
                % - Data Structure
                % - Results Structure
                %

                % Revision History
                %
                % Versions now tracked on GitHub https://github.com/Mensen/ept_TFCE-matlab
                %
                % Version 2.2
                % 12.09.2013
                % - Calls the new TFCE mex file which simultaneously looks for negative
                % values (should be significantly faster then two separate calls)
                %
                % Version 2.1
                % 11.09.2013
                % - Included if statements for 2D and 3D TFCE differentiation
                %
                % 22.01.2013
                % Corrected bug in correlation analysis (check for type 'r' instead of 'c')
                % Eliminated the use of the shuffle function for correlation
                % 
                % 06.11.2012
                % Included a one-sample test and some error handling
                %
                % 14.10.2012
                % Includes Info/Data/Results Structure for better overview of test
                % Includes the type 'c' for correlational analysis
                % 23.12.2010
                % - Performs two TFCE calculations, positive and negative values and then
                %   puts them back together
                %       - To prevent non existent clusters forming when taking abs values
                %
                % 03.01.2011
                % - Calculated significance power and displays graph
                %       - Significance power is the inverse of the p-value squared...
                % 
                % 10.03.2011
                % - Now a function that can be called from the command line
                %
                % 24.05.2011
                % - Added random stream dependence on the clock
                % 
                % 12.07.2011
                % - Adapted to loading the two summary files to be compared
                % - Uses new modified version of the channel neighbours algorithm
                %
                % 16.12.2011
                % - Used -log of p-values to plot (more interpretable than power values)
                % - Assigns defaults rather than prompting user for information
                %  







                % assignin('base', 'varargin', varargin); 
                %% Set Defaults

                 E_H        = [0.66 2]; % default parameters of E and H
                 nPerm      = 5000; % default number of permutations
                 rSample    = 250; % default assumed sampling rate (not used in calculation but only for result viewer)
                 saveName   = ['ept_Results_', date, '.mat'];
                 type       = 'i'; % i = independent sample T-Test; d = dependent (paired) sample T-Test; c = Pearson's correlation
                 plots      = 0; % change to '1' to show significance plots after calculation
                 flag_ft    = 0; 
                 flag_tfce  = 1;
                 flag_save  = 1;
                 ChN = [];

                % Set random stream depending on the clock value (unpredictable).
                myRand = RandStream('mt19937ar','Seed',sum(100*clock));
                RandStream.setGlobalStream(myRand);

                % Process Secondary Arguments
                if nargin > 2
                  if (round(nargin/2) == nargin/2)
                    error('Even number of input arguments??')
                  end
                  for i = 1:2:length(varargin)
                    Param = varargin{i};
                    Value = varargin{i+1};
                    if ~ischar(Param)
                      error('Flag arguments must be strings')
                    end
                    Param = lower(Param);

                    switch Param
                        case 'e_h'
                            E_H         = Value;
                        case 'nperm'
                            nPerm       = Value;
                        case 'rsample'
                            rSample     = Value;
                        case 'fsample'
                            fSample     = Value;
                        case 'savename'
                            saveName    = Value;
                        case 'type'
                            type        = lower(Value);
                        case 'flag_ft'
                            flag_ft     = Value;
                        case 'flag_tfce'
                            flag_tfce   = Value;
                        case 'flag_save'
                            flag_save   = Value;
                        case 'plots'
                            plots       = Value;
                        case 'chn'
                            ChN         = Value;
                        otherwise
                            display (['Unknown parameter setting: ' Param])
                    end
                  end
                end

                %% Process Arguments
                if nargin < 1;

                    prompt = {'Number of Permutations', 'Sampling Rate', '(i)n(d)ependent / (c)orrelation?', 'Results Name', 'Plots? 0=No, 1=Yes'};
                    dlg_title = 'Define Parameters';
                    def = {num2str(nPerm), num2str(rSample), type, saveName, num2str(plots)};
                    noGr = inputdlg(prompt, dlg_title,1,def);

                    if isempty(noGr)
                        error('Analysis cancelled')
                    end

                    nPerm      = str2double(noGr{1});  % number of permutations
                    rSample    = str2double(noGr{2});  % sampling rate
                    type       = noGr{3};
                    saveName   = noGr{4};
                    plots      = str2double(noGr{5});

                    if type == 'i' || type == 'd'
                        [DataFile{1}, DataPath{1}] = uigetfile('', 'Please Select the first EEG Summary File', 'MultiSelect', 'off');
                        [DataFile{2}, DataPath{2}] = uigetfile('', 'Please Select the second EEG Summary File', 'MultiSelect', 'off');
                    elseif type == 'c'
                        [DataFile{1}, DataPath{1}] = uigetfile('', 'Please Select the EEG Summary File', 'MultiSelect', 'off');
                        [DataFile{2}, DataPath{2}] = uigetfile('', 'Please Select the Behavioural Summary File', 'MultiSelect', 'off');
                    elseif type == 'o'
                        [DataFile{1}, DataPath{1}] = uigetfile('', 'Please Select the EEG Summary File', 'MultiSelect', 'off');
                    end

                    if DataFile{1} == 0
                        fprintf(1, 'Returning: No Data File Selected... \n');
                        return;
                    end

                    FullFileName1 = strcat(DataPath{1}, DataFile{1});
                    Data{1,1}     = load (FullFileName1);
                    Data{1,1}     = Data{1}.Summary;

                    % Error handling for non-onesample tests...
                    if type ~= 'o'
                        FullFileName2 = strcat(DataPath{2}, DataFile{2});
                        Data{2,1}     = load (FullFileName2);
                    else
                        Data{2,1}     = zeros(size(Data{1}));
                    end

                    % Error handling for special behavioural data exception
                    if type == 'i' || type == 'd'
                        Data{2,1}	= Data{2}.Summary;
                        aData = [Data{1};Data{2}];
                    elseif type == 'c'
                        Data{2}     = Data{2}.Behavioural;
                    end

                    % Load the electrode file...
                    [ElecFile, ElecPath] = uigetfile('', 'Please Select the Electrodes File');
                    FullElecName = strcat(ElecPath, ElecFile);
                    load (FullElecName);

                elseif nargin > 2;
                    % check if file or variable entered as data
                    if isa(DataFile1, 'char')
                        % Not a lot of error checking done here yet.
                        DataFile{1} = DataFile1;
                        DataFile{2} = DataFile2;

                        Data{1}       = load (DataFile1, 'Summary');
                        Data{1}       = Data{1}.Summary;
                        Data{2}       = load (DataFile2);
                        
                        if type == 'i' || type == 'd'
                            Data{2} 	= Data{2}.Summary;
                            aData = [Data{1};Data{2}];
                        elseif type == 'c'
                            Data{2}     = Data{2}.Behavioural;
                        end

                    elseif isa(DataFile1, 'double')
                        DataFile{1, 2} = '';
                        Data{1} = DataFile1;
                        Data{2} = DataFile2;
                        if type == 'i' || type == 'd'
                            aData = [Data{1}; Data{2}];
                        end

                    elseif isa(DataFile1, 'single')
                        fprintf(1, 'warning: single data converted to double for analysis');
                        DataFile{1, 2} = '';
                        Data{1} = double(DataFile1);
                        Data{2} = double(DataFile2);
                        if type == 'i' || type == 'd'
                            aData = [Data{1}; Data{2}];
                        end 

                    else
                        fprintf(1, 'error: could not recognise data format (possibly single) \n');
                    end

                    if isa(ElecFile, 'char')
                        e_loc         = load (ElecFile);
                        e_loc         = e_loc.e_loc;
                    else
                        e_loc = ElecFile;
                    end
                    
                end

                nA   = size(Data{1},1);
                nB   = size(Data{2},1);
                nCh  = size(Data{1},2);
                nS   = size(Data{1},3);

                % For Frequency-Time Data...
                if ndims(Data{1})==4;
                    nS = size(Data{1},4);
                    nF = size(Data{1},3);
                else
                    nF = [];
                end

                %% -- Error Checking -- %%
                if ~isequal(nB, nA) && type == 'r'
                    error ('Must have the same number of participants as behavioural points...')
                end

                if ~isequal(nB, nA) && type == 'd'
                    error ('Must have the same number of participants for paired comparisons...')
                end

                % Check Location File for Number of Channels
                if ~flag_ft
                    if ~isequal(nCh, length(e_loc))
                        error ('Number of channels in data does not equal that of locations file')
                    end
                end

                %% Create Summary File

                % Summary = [Data{1};  Data{2}];
                tic; % Start the timer for the entire analysis

                %% Calculate the channels neighbours... using the modified version ChN2

                if ~flag_ft && isempty(ChN)
                    display('Calculating Channel Neighbours...')
                    ChN = ept_ChN2(e_loc);
                    display('Done')
                end

                %% Create all variables in loop at their maximum size to increase performance

                maxTFCE = zeros(nPerm,1);

                %% Calculate the actual T values of all data

                display('Calculating Actual Differences...')

                % Calculate different T-values for independent and dependent groups
                if type == 'i'

                    T_Obs = (mean(Data{1})-mean(Data{2}))./sqrt(var(Data{1})/nA+var(Data{2})/nB);
                    T_Obs = squeeze(T_Obs);

                elseif type == 'd' || type == 'o'

                    D    = Data{1}-Data{2};

                    T_Obs = mean(D,1)./(std(D)./sqrt(nA));
                    T_Obs = squeeze(T_Obs);

                elseif type == 'c'
                    % Define "condition" function for correlation calculation
                    condition = @(x) (x-mean(x))./std(x);

                    % calculate the size of the data
                    data_size = size(Data{1}); data_size(1) = [];

                    % repmat the correlation data to match the data size
                    behavior_repeated = repmat(Data{2}, [1, data_size]);

                    % calculate the correlation coefficient (r) but call it T_Obs for consistency
                    temp_observed = [reshape(condition(Data{1}), nA, [])' ...
                        * reshape(condition(behavior_repeated), nA, [])] ...
                        / reshape(sum(condition(Data{1}).^2), 1, []);
                    T_Obs = reshape(temp_observed, data_size);

                end

                % check for non-zero T-values (will crash TFCE)
                if max(abs(T_Obs(:))) < 0.00001
                    error('T-values were all 0')
                end

                % TFCE transformation...
                if flag_tfce
                    if ismatrix(T_Obs);
                        T_Obs=double(T_Obs);%agragado magico
                        if ~flag_ft
                            TFCE_Obs = ept_mex_TFCE2D(T_Obs, ChN, E_H);
                        else
                            % if data is not in channel neighbourhood
                            % artificially create 3rd dimension
                            T_Obs = repmat(T_Obs, [1, 1, 2]);
                            deltaT = max(abs(T_Obs(:)))/50;
                            TFCE_Obs = ept_mex_TFCE(T_Obs, deltaT);
                            % remove extra dimension
                            T_Obs = T_Obs(:, :, 1);
                            TFCE_Obs = TFCE_Obs(:, :, 1);
                        end
                    end

                    if ndims(T_Obs) == 3;
                        TFCE_Obs = ept_mex_TFCE3D(T_Obs, ChN, E_H);
                    end
                else
                    TFCE_Obs = T_Obs;
                end
                display('Done')

                %% Calculating the T value and TFCE enhancement of each different permutation

                display('Calculating Permutations...')

                    for i   = 1:nPerm           % Look into parfor for parallel computing

                            if type == 'i' %two sample independent T-test
                                r_perm      = randperm(nA+nB); % Consider using Shuffle mex here (50-85% faster)...

                                if ismatrix(T_Obs);
                                    nData       = aData(r_perm,:,:); 
                                    sData{1}    = nData(1:nA,:,:); sData{2}= nData((nA+1):(nA+nB),:,:);
                                else
                                    nData       = aData(r_perm,:,:,:); 
                                    sData{1}    = nData(1:nA,:,:,:); sData{2}= nData((nA+1):(nA+nB),:,:,:);
                                end 

                                T_Perm = (mean(sData{1})-mean(sData{2}))./sqrt(var(sData{1})/nA+var(sData{2})/nB);
                                T_Perm = squeeze(T_Perm);   

                            elseif type == 'd' || type == 'o' %one-sample T-test
                                Signs =[-1,1];
                                SignSwitch = randsample(Signs,size(D,1),'true')';
                                if ismatrix(T_Obs);
                                   SignSwitch = repmat(SignSwitch, [1 nCh nS]);
                                else
                                   SignSwitch = repmat(SignSwitch, [1 nCh nF nS]);
                                end 

                                nData = SignSwitch.*D;

                                T_Perm = mean(nData,1)./(std(nData)./sqrt(size(nData,1)));
                                T_Perm = squeeze(T_Perm);


                            elseif type == 'c' %correlation analysis

                                % repmat the correlation data to match the data size
                                behavior_perm = Data{2}(randperm(nB));
                                behavior_repeated = repmat(behavior_perm, [1, data_size]);

                                % calculate the correlation coefficient (r) but call it T_Perm for consistency
                                temp_observed = [reshape(condition(Data{1}), nA, [])' ...
                                    * reshape(condition(behavior_repeated), nA, [])] ...
                                    / reshape(sum(condition(Data{1}).^2), 1, []);
                                T_Perm = reshape(temp_observed, data_size);


                            else
                                error('Unrecognised analysis-type; see help file for valid inputs')
                            end

                        % TFCE transformation...
                        if flag_tfce
                            if ismatrix(T_Perm)
                                T_Perm = double(T_Perm);
                                if ~flag_ft
                                    TFCE_Perm = ept_mex_TFCE2D(T_Perm, ChN, E_H);
                                else
                                    % if data is not in channel neighbourhood
                                    % artificially create 3rd dimension
                                    T_Perm = repmat(T_Perm, [1, 1, 2]);
                                    deltaT = max(abs(T_Perm(:)))/50;
                                    TFCE_Perm = ept_mex_TFCE(T_Perm, deltaT);
                                    % remove extra dimension
                                    T_Perm = T_Perm(:, :, 1);
                                    TFCE_Perm = TFCE_Perm(:, :, 1);
                                end
                            end

                            if ndims(T_Perm) == 3;
                                TFCE_Perm = ept_mex_TFCE3D(T_Perm, ChN, E_H);
                            end

                        else
                            TFCE_Perm = T_Perm;
                        end

                        maxTFCE(i) = max(abs(TFCE_Perm(:)));       % stores the maximum absolute value

                        progressbar(i/nPerm); %Progress Bar

                    end

                display('Done')

                %% Calculating the p value from the permutation distribution

                display('Calculating P-Values and Saving...')

                % add observed maximum
                edges = [maxTFCE;max(abs(TFCE_Obs(:)))];

                [~,bin]     = histc(abs(TFCE_Obs),sort(edges));
                P_Values    = 1-bin./(nPerm+2);

                % plot the "significance power"...

                if plots == 1;
                    if ismatrix(T_Perm);
                        figure
                        plot(sum(-log(P_Values)));
                        title (saveName)
                    end
                end

                % Save test information in single structure
                c = clock;
                Info.Comments = ['TFCE analysis conducted at ', num2str(c(4)), ':', num2str(c(5)), ' on ', date];

                Info.Parameters.E_H         = E_H;
                Info.Parameters.nPerm       = nPerm;
                Info.Parameters.rSample     = rSample;
                Info.Parameters.type        = type;
                Info.Parameters.nChannels   = nCh;
                Info.Parameters.nSamples    = nS; % Not sure if actually used...
                Info.Parameters.GroupSizes  = [nA, nB];

                if exist('fSample', 'var')
                    Info.Parameters.fSample = fSample;
                end

                Info.Electrodes.e_loc       = e_loc;
                Info.Electrodes.ChannelNeighbours = ChN;

                Info.DataFiles              = DataFile;

                Results.Obs                 = T_Obs;
                Results.TFCE_Obs            = TFCE_Obs;
                Results.maxTFCE             = sort(maxTFCE);
                Results.P_Values            = P_Values;

                % save the file in the current directory
                if flag_save
                    save (saveName, 'Info', 'Data', 'Results', '-mat')
                end

                %%
                display('All done!')
                toc

                [min_P, idx] = min(Results.P_Values(:));
                [Ch, S]      = ind2sub(size(Results.P_Values),idx);
                max_Obs      = Results.Obs(idx);

                display(['Peak significance found at channel ', num2str(Ch), ' at sample ', num2str(S), ': T(', num2str(size(Data{1},1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);


        end
        function [distractors,targets,easyHard]         = fun_conditionHistogramInd(EEG,fixsSuj,minfixdur) 
            %event = EEG.event;
            fixs = fixsSuj;
            %keyboard
            %urevents = arrayfun(@(y) y.eventurevent{cellfun(@(x) x==0, y.eventlatency)},EEG.epoch); %busco eventos de la fijacion principal
            infix = ismember([fixs.urevent],[fixs.urevent]); % indices de las fijaciones que corresponden al fixs con urevents del epoqueado por fixs



            %index to retain
            distractors = [];
            targets     = [];
            easyHard    = [];
            
            tmax = minfixdur ;
            indfixnum =  find(infix);
            for i = 1:sum(infix)
                    in = indfixnum(i);
                    %urin   = find([urevents]== fixs(in).urevent);
                    if (  (fixs(in).refix==0 || fixs(in).refix==1)  && strcmp(fixs(in).trial_cattype,'F') &&... strcmp(fixs(in).trial_type,'VS') &&... strcmp(fixs(in).trial_cattype,'I') &&... 
                            fixs(in).dur_added>tmax && ...
                            strcmp(fixs(in).distractor_cat,'faces')) ...&& strcmp(fixs(in).trial_cat,'faces'))
                            %indexs to retain
                            distractors = [distractors fixs(in).dur_added];
                    end

                    if ( (fixs(in).refix==0 || fixs(in).refix==1) &&...  strcmp(fixs(in).trial_type,'VS') &&... strcmp(fixs(in).trial_cattype,'I') &&...
                            fixs(in).dur_added>tmax && ...
                            strcmp(fixs(in).distractor_cat,'objects'))... && strcmp(fixs(in).trial_cat,'faces'))
                            %indexs to retain

                            distractors = [distractors fixs(in).dur_added];
                    end
                    if ( (fixs(in).refix==0 || fixs(in).refix==1) && strcmp(fixs(in).trial_type,'VS') &&...
                            fixs(in).dur_added>tmax && fixs(in).istarget==1 && fixs(in).trial_correct == 1 &&...
                            strcmp(fixs(in).trial_resp,'P') && strcmp(fixs(in).trial_cattype,'F'))
                            %indexs to retain

                            targets = [targets fixs(in).dur_added];
                   
                    end
                    if ( (fixs(in).refix==0 || fixs(in).refix==1) && strcmp(fixs(in).trial_type,'VS') && ...
                            fixs(in).dur_added>tmax && fixs(in).istarget==1 && fixs(in).trial_correct == 1 &&...
                            strcmp(fixs(in).trial_resp,'P') && strcmp(fixs(in).trial_cat,'faces') && strcmp(fixs(in).trial_cattype,'I')  )
                            %indexs to retain

                            targets = [targets fixs(in).dur_added];
                    end
                     if (  (fixs(in).refix==0 || fixs(in).refix==1) && strcmp(fixs(in).trial_type,'VS') && ...&& strcmp(fixs(in).trial_cattype,'O') &&... 
                            fixs(in).dur_added>tmax && ...
                            strcmp(fixs(in).distractor_cat,'objects') && strcmp(fixs(in).trial_cat,'faces'))
                            %indexs to retain
                            easyHard = [easyHard fixs(in).dur_added];
                    
                    end
                    if ( (fixs(in).refix==0 || fixs(in).refix==1) &&  strcmp(fixs(in).trial_type,'VS') &&... strcmp(fixs(in).trial_cattype,'F') &&...
                            fixs(in).dur_added>tmax && ...
                            strcmp(fixs(in).distractor_cat,'objects') && strcmp(fixs(in).trial_cat,'objects'))
                            %indexs to retain

                            easyHard = [easyHard fixs(in).dur_added];
                    end                              
    
            end
        end
        function EEG                                    = renameEventsUnfold(EEG, fixsSuj,minfixdur,istarget)    
            event = EEG.event;
            fixs = fixsSuj;
            
            %EEG = pop_rmbase( EEG, [-199.2188 0]);%baseline to 0 secs 27/02/2020 
            
            %keyboard
            urevents = [event(ismember({event.type},'fixation')).urevent]; %urevents from EEG corresponding to fixations
            infix = ismember([fixs.urevent],urevents); % mask with indices in Fixs that corresponds to urevents in EEG im gonna iterate over them.
            %this is to exclude events that could no belong to the EEG
            %already epoched by trials
            
            marc = {'13',  'R_blink' ,  'R_fixation' , 'R_saccade' ,   'bad_ET'   , 'boundary', 'saccade','L_blink','L_fixation','L_saccade'}
            EEG.event(ismember({EEG.event.type}, marc)) = []
                
             
            tmax = minfixdur ;
            indfixnum =  find(infix);
            for i = 1:sum(infix)
                    in = indfixnum(i);
                    urin   = find([EEG.event.urevent]== fixs(in).urevent); %tell me which index of event corresponds to this urevent from fix
                    %keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  if (  (fixs(in).refix==0 || fixs(in).refix==1) && strcmp(fixs(in).trial_type,'VS') && ...&& strcmp(fixs(in).trial_cattype,'O') &&... 
%                             fixs(in).dur_added>tmax && ...
%                             strcmp(fixs(in).distractor_cat,'objects') && strcmp(fixs(in).trial_cat,'faces'))
%                             %indexs to retain
%                             NTO_easy = [NTO_easy urin];
%                     
%                     end
%                     if ( (fixs(in).refix==0 || fixs(in).refix==1) &&  strcmp(fixs(in).trial_type,'VS') &&... strcmp(fixs(in).trial_cattype,'F') &&...
%                             fixs(in).dur_added>tmax && ...
%                             strcmp(fixs(in).distractor_cat,'objects') && strcmp(fixs(in).trial_cat,'objects'))
%                             %indexs to retain
% 
%                             NTO_hard = [NTO_hard urin];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                   
                    if    ( (fixs(in).refix==0 || fixs(in).refix==1) && strcmp(fixs(in).trial_type,'VS') &&...%agrego VS t para analizar easy
                                fixs(in).dur_added>tmax && fixs(in).istarget == istarget  && fixs(in).dur_added< 1)%...agregue corte a 1seg
                                %&& fixs(in).trial_correct==1)
                         EEG.event(urin).istarget = fixs(in).istarget;
                         EEG.event(urin).fixationRank = fixs(in).n;
                        
                        if      ( strcmp(fixs(in).distractor_cat,'faces') )
                                EEG.event(urin).stimulusType = 1;                   %'face';
                                EEG.event(urin).stimulusDur = fixs(in).dur_added;
                             
                                %EEG.event(urin).(['Subj' num2str(subj)]) = 1;
                                EEG.event(urin).type = 'erase';
                                %agrego para easy; borra las caras                             
                        elseif   ( strcmp(fixs(in).distractor_cat,'objects'))
                                EEG.event(urin).stimulusType = 0;%'object';
                                EEG.event(urin).stimulusDur = fixs(in).dur_added;
                            
                                %EEG.event(urin).(['Subj' num2str(subj)]) = 1;                      
                        end
                        
                        if   strcmp(fixs(in).trial_cat,'faces')
                                EEG.event(urin).trialType    = 'face';
                        % elseif   strcmp(fixs(in).trial_cattype,'I') 
                         %       EEG.event(urin).type = 'erase';
                        %else
                        %     EEG.event(urin).type = 'erase';%agrego para easy
                        end
                        if  strcmp(fixs(in).trial_cat,'objects')
                                EEG.event(urin).trialType    = 'object';
                        end
                     
                    else
                         EEG.event(urin).type = 'erase';
                       

                 
                    end
                    %if strcmp(fixs(in).trial_cattype,'I')
                     %   EEG.event(urin).type = 'erase';
                    %end
                    if fixs(in).isdistractor == 0
                        EEG.event(urin).type = 'erase';
                    end
            end
            EEG.event(~ismember(urevents,[fixs.urevent])) =[];
            EEG.event(ismember({EEG.event.type},'erase')) = [];
        end 
        function EEG                                    = renameEventsUnfold2021(EEG, TSuj,hipothesis,minfixdur)    
            event   = EEG.event;
            fixs    = TSuj;
            
            %EEG = pop_rmbase( EEG, [-199.2188 0]);%baseline to 0 secs 27/02/2020 

            urevents= [event(ismember({event.type},'fixation')).urevent]; %urevents from EEG corresponding to fixations
            infix   = ismember([fixs.urevent],urevents); % mask with indices in Fixs that corresponds to urevents in EEG im gonna iterate over them.
            %this is to exclude events that might belong to the EEG
            %already epoched by trials
            marc = {'13',  'R_blink' ,  'R_fixation' , 'R_saccade' ,  'saccade','L_blink','L_fixation','L_saccade'};
            EEG.event(ismember({EEG.event.type}, marc)) = [];

            %n=0;
            tmax = minfixdur ;
            indfixnum =  find(infix);
            for in = indfixnum'

                    urin   = find([EEG.event.urevent] == fixs.urevent(in)); %tells which index of event corresponds to this urevent from fix
                    %keyboard
                    if    ( (fixs.refix(in) ==0 || fixs.refix(in) ==1) && fixs.dur_added(in) >tmax && ...
                             fixs.dur_added(in)<1.5)%...agregue corte a 1seg % mdf 22/6/21 borre:  &&  fixs.isdistractor(in) == 1 
                         %  && ~isnan(fixs.(hipothesis)(in))
                       
                        % n = n +1            
                        % keyboard
                        EEG.event(urin).fixationRank    = fixs.rank(in);            % mdf 22/6/21 raw rank
                        EEG.event(urin).StandFixRank    = fixs.Nrank(in);           % mdf 22/6/21 normalized rank [(rank-min_rank)/(max_rank-min_rank)], min_rank and max_rank per trial
                        EEG.event(urin).GlobalNRank     = fixs.N_subj_rank(in);     % mdf 1/4/22 normalized rank [(rank-min_rank)/(max_rank-min_rank)], min_rank and max_rank per subject
                        EEG.event(urin).NCRank          = fixs.NC_rank(in);         % mdf 1/4/22 StandFixRank centered [StandFixRank-0.5]
                        EEG.event(urin).GlobalNCRank    = fixs.NC_subj_rank(in);    % mdf 1/4/22 GlobalNRank centered [GlobalNRank-0.5]
                        EEG.event(urin).isEX            = fixs.isEX(in);
                        EEG.event(urin).catRank         = ~fixs.first_rank(in);     % 1 for fixationRank less than 5 % mdf 26/7/21 
                        EEG.event(urin).(hipothesis)    = fixs.(hipothesis)(in);  
                        EEG.event(urin).stimulusDur     = fixs.dur_added(in);
                        % EEG.event(urin).saccade_amp = fixs.saccade_amplitud(in); 

                        %DAC mod oct 2021
                        if fixs.isdistractor(in)
                            EEG.event(urin).faces       = strcmp(fixs.distractor_cat(in),'faces');
                        elseif  fixs.istarget(in)
                            EEG.event(urin).faces       = strcmp(fixs.trial_cat(in),'faces');
                        end
                    else
                         EEG.event(urin).type = 'erase';              
                    end
            end

            for ur = [EEG.event.urevent] 
                if ~ismember(ur,[fixs.urevent]) && ~strcmp({EEG.event(find([EEG.event.urevent]==ur)).type},'bad_ET')
                    EEG.event(find([EEG.event.urevent]==ur)) =[];
                end
            end
            EEG.event(ismember({EEG.event.type},'erase')) = [];
        end 
        function EEG                                    = renameEventsUnfoldSaccade(EEG, TSuj,hipothesis,minfixdur)    
            event   = EEG.event;
            fixs    = TSuj;
            
            %EEG = pop_rmbase( EEG, [-199.2188 0]);%baseline to 0 secs 27/02/2020 

            urevents= [event(ismember({event.type},'fixation')).urevent]; %urevents from EEG corresponding to fixations
            infix   = ismember([fixs.urevent],urevents); % mask with indices in Fixs that corresponds to urevents in EEG im gonna iterate over them.
            %this is to exclude events that might belong to the EEG
            %already epoched by trials
%             marc = {'13',  'R_blink' ,  'R_fixation' , 'R_saccade' ,  'saccade','L_blink','L_fixation','L_saccade'};
            marc = {'13','R_blink','R_fixation','R_saccade','L_blink','L_fixation','L_saccade'};
            EEG.event(ismember({EEG.event.type}, marc)) = [];

            %n=0;
            tmax = minfixdur ;
            indfixnum =  find(infix);
            for in = indfixnum'

                    urin   = find([EEG.event.urevent] == fixs.urevent(in)); %tells which index of event corresponds to this urevent from fix
                    %keyboard
                    if    ( (fixs.refix(in) ==0 || fixs.refix(in) ==1) && fixs.dur_added(in) >tmax && ...
                             fixs.dur_added(in)<1.5)%...agregue corte a 1seg % mdf 22/6/21 borre:  &&  fixs.isdistractor(in) == 1 
                         %  && ~isnan(fixs.(hipothesis)(in))
                       
                        % n = n +1            
                        % keyboard
                        EEG.event(urin).fixationRank    = fixs.rank(in);            % mdf 22/6/21 raw rank
                        EEG.event(urin).StandFixRank    = fixs.Nrank(in);           % mdf 22/6/21 normalized rank [(rank-min_rank)/(max_rank-min_rank)], min_rank and max_rank per trial
                        EEG.event(urin).GlobalNRank     = fixs.N_subj_rank(in);     % mdf 1/4/22 normalized rank [(rank-min_rank)/(max_rank-min_rank)], min_rank and max_rank per subject
                        EEG.event(urin).NCRank          = fixs.NC_rank(in);         % mdf 1/4/22 StandFixRank centered [StandFixRank-0.5]
                        EEG.event(urin).GlobalNCRank    = fixs.NC_subj_rank(in);    % mdf 1/4/22 GlobalNRank centered [GlobalNRank-0.5]
                        EEG.event(urin).isEX            = fixs.isEX(in);
                        EEG.event(urin).catRank         = ~fixs.first_rank(in);     % 1 for fixationRank less than 5 % mdf 26/7/21 
                        if ~strcmp(hipothesis,'faces') % (JK, 2023-04-20): This is fix below
                            EEG.event(urin).(hipothesis)    = fixs.(hipothesis)(in);  
                        end
                        EEG.event(urin).stimulusDur     = fixs.dur_added(in);
                        % EEG.event(urin).saccade_amp = fixs.saccade_amplitud(in); 

                        %DAC mod oct 2021
                        if fixs.isdistractor(in)
                            EEG.event(urin).faces       = strcmp(fixs.distractor_cat(in),'faces');
                        elseif  fixs.istarget(in)
                            EEG.event(urin).faces       = strcmp(fixs.trial_cat(in),'faces');
                        end
                    else
                         EEG.event(urin).type = 'erase';              
                    end
            end

            for ur = [EEG.event.urevent] 
                if ~ismember(ur,[fixs.urevent]) && ...
                        ~strcmp({EEG.event(find([EEG.event.urevent]==ur)).type},'bad_ET')  && ...
                        ~strcmp({EEG.event(find([EEG.event.urevent]==ur)).type},'saccade')
                    EEG.event(find([EEG.event.urevent]==ur)) =[];
                end
            end
            EEG.event(ismember({EEG.event.type},'erase')) = [];
        end 
        function plotUnfoldErpimageCheck(EEG,channel,hipothesis,alignto,win)

            chi = find(ismember({EEG.chanlocs.labels},channel));

            cfg = [];
            %cfg.title = 'nombre';
            cfg.erp = 'on' ;
            cfg.split_by = hipothesis;
            cfg.channel = chi;
            cfg.figure = 1;
            cfg.caxis = [-5,5];
            cfg.alignto = alignto;
            cfg.win = win;
            cfg = [fieldnames(cfg),struct2cell(cfg)].';
            cfg(:)

            % uf_erpimage_mod(EEG,'type','raw','sort_by','duration',cfg{:});
            % %%
            % uf_erpimage_mod(EEG,'type','modelled','sort_by','NTF_NTO',cfg{:});
            uf_erpimage_mod(EEG,'type','raw','sort_by','duration',cfg{:});

            uf_erpimage_mod(EEG,'type','modelled','sort_by','duration','addResiduals', 0,cfg{:});

            uf_erpimage_mod(EEG,'type','modelled','sort_by','duration','addResiduals', 1,cfg{:});
       end
        function EEG                                    = filterband(inFile, locutoff, hicutoff,minphase)
            EEG = pop_loadset(inFile);
            EEG = eeg_checkset(EEG);
            order = 3*fix(EEG.srate/locutoff);
            EEG = pop_eegfiltnew(EEG, locutoff, hicutoff,[],false,[],[],minphase,[]);
            for ch=1:size(EEG.data,1)
                EEG.data(ch,:,:) = abs(hilbert(EEG.data(ch,:,:)));
            end
        end

    end
end

%% Functions Definition
function name = printLoadingFile(fullPath)
    % check version name = split(EDF_file,'/');
    name = strsplit(fullPath,'/');
    name =name{end}
end
function [ChN] = ept_ChN2(eLoc, displayNet)
% Searches for neighbouring channels using triangulation and reports back each channels neighbours

% This file is part of the program ept_TFCE.
% ept_TFCE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% ept_TFCE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with ept_TFCE.  If not, see <http://www.gnu.org/licenses/>.

% 11.12.2013
% Fixed error in script for small electrode files where channel list was
% not cleared properly in each iteration of finding neighbours...

if nargin < 2
    displayNet = 0; %set to 1 if you want net to be displayed in a figure.
end

nCh = length(eLoc);

indices = 1:nCh;

x = [eLoc(indices).X ]'; % Gets the X,Y, and Z values from loc_file
y = [eLoc(indices).Y ]';
z = [eLoc(indices).Z ]';

vertices = [x,y,z];

%% [X,Y] = bst_project_2d(vertices(:,1), vertices(:,2), vertices(:,3));

        z2 = z - max(z);

	%% [TH,PHI,R] = cart2sph(x, y, z) // [az,elev,r] = cart2sph(x,y,z)

        hypotxy = hypot(x,y);
        R = hypot(hypotxy,z2);
        PHI = atan2(z2,hypotxy);
        TH = atan2(y,x);

        % Remove the too small values for PHI
        PHI(PHI < 0.001) = 0.001;

        % Flat projection
        R2 = R ./ cos(PHI) .^ .2;

    %% [X,Y] = pol2cart(TH,R2) // [x,y,z] = pol2cart(th,r,z)

        X = R2.*cos(TH);
        Y = R2.*sin(TH);

%% bfs_center = bst_bfs(vertices)' // [ HeadCenter, Radius ] = bst_bfs( Vertices )

    mass = mean(vertices);
    diffvert = bsxfun(@minus, vertices, mass); % originally uses bst_bsxfun but only for older Matlab versions
    R0 = mean(sqrt(sum(diffvert.^2, 2)));
    % Optimization
    vec0 = [mass,R0];
    minn = fminsearch(@dist_sph, vec0, [], vertices);
    HeadCenter = minn(1:end-1); % 3x1

    
coordC = bsxfun(@minus, vertices, HeadCenter);
coordC = bsxfun(@rdivide, coordC, sqrt(sum(coordC.^2,2)));
coordC = bsxfun(@rdivide, coordC, sqrt(sum(coordC.^2,2)));
% Tesselation of the sensor array
faces  = convhulln(coordC);


%% Remove unnecessary triangles...

% Get border of the representation
border = convhull(X,Y);
%plot(X(border),Y(border),'r-',X,Y,'b+')

% Keep faces inside the border
iInside = ~(ismember(faces(:,1),border) & ismember(faces(:,2),border)& ismember(faces(:,3),border));
faces   = faces(iInside, :);

    my_norm = @(v)sqrt(sum(v .^ 2, 2)); % creates an object
    % Get coordinates of vertices for each face
    vertFacesX = reshape(vertices(reshape(faces,1,[]), 1), size(faces));
    vertFacesY = reshape(vertices(reshape(faces,1,[]), 2), size(faces));
    vertFacesZ = reshape(vertices(reshape(faces,1,[]), 3), size(faces));
    % For each face : compute triangle perimeter
    triSides = [my_norm([vertFacesX(:,1)-vertFacesX(:,2), vertFacesY(:,1)-vertFacesY(:,2), vertFacesZ(:,1)-vertFacesZ(:,2)]), ...
                my_norm([vertFacesX(:,1)-vertFacesX(:,3), vertFacesY(:,1)-vertFacesY(:,3), vertFacesZ(:,1)-vertFacesZ(:,3)]), ...
                my_norm([vertFacesX(:,2)-vertFacesX(:,3), vertFacesY(:,2)-vertFacesY(:,3), vertFacesZ(:,2)-vertFacesZ(:,3)])];
    triPerimeter = sum(triSides, 2);
    % Threshold values
    thresholdPerim = mean(triPerimeter) + 3 * std(triPerimeter);
    % Apply threshold
    faces(triPerimeter > thresholdPerim, :) = [];

    
%% Display Net
if displayNet == 1;

    figure( 'Color',       'w'     ,...
                'Position',    [50,50, 500, 500]  );

        axes  ( 'Color',      'w' );   

        FaceColor = [.5 .5 .5];
        EdgeColor = [0 0 0];
        FaceAlpha = .9;
        LineWidth = 1;

        hNet = patch('Vertices',        vertices, ...
                     'Faces',           faces, ...
                     'FaceVertexCData', repmat([1 1 1], [length(vertices), 1]), ...
                     'Marker',          'o', ...
                     'LineWidth',       LineWidth, ...
                     'FaceColor',       FaceColor, ...
                     'FaceAlpha',       FaceAlpha, ...
                     'EdgeColor',       EdgeColor, ...
                     'EdgeAlpha',       1, ...
                     'MarkerEdgeColor', [0 0 0], ...
                     'MarkerFaceColor', 'flat', ...
                     'MarkerSize',      12, ...
                     'BackfaceLighting', 'lit', ...
                     'Tag',             'SensorsPatch');
        material([ 0.5 0.50 0.20 1.00 0.5 ])
        lighting phong
    % Set Constant View Angle
    view(48,15);
    axis equal
    cameratoolbar('Show')
    rotate3d on
    axis off
       
end

% loop through all the channels
output = cell(nCh,1);
for n = 1 : nCh
   
    [r1, ~] = ind2sub(size(faces), find(faces == n));
    output{n} = unique(faces(r1, :)')';
    
end

nz=max(cellfun(@numel,output));
ChN=cell2mat(cellfun(@(a) [a,zeros(1,nz-numel(a))],output,'uni',false));

end
function d = dist_sph(vec,sensloc)
R = vec(end);
center = vec(1:end-1);
% Average distance between the center if mass and the electrodes
diffvert = bsxfun(@minus, sensloc, center);
d = mean(abs(sqrt(sum(diffvert.^2,2)) - R));
end
function progressbar(varargin)
% Description:
%   progressbar() provides an indication of the progress of some task using
% graphics and text. Calling progressbar repeatedly will update the figure and
% automatically estimate the amount of time remaining.
%   This implementation of progressbar is intended to be extremely simple to use
% while providing a high quality user experience.
%
% Features:
%   - Can add progressbar to existing m-files with a single line of code.
%   - Supports multiple bars in one figure to show progress of nested loops.
%   - Optional labels on bars.
%   - Figure closes automatically when task is complete.
%   - Only one figure can exist so old figures don't clutter the desktop.
%   - Remaining time estimate is accurate even if the figure gets closed.
%   - Minimal execution time. Won't slow down code.
%   - Randomized color. When a programmer gets bored...
%
% Example Function Calls For Single Bar Usage:
%   progressbar               % Initialize/reset
%   progressbar(0)            % Initialize/reset
%   progressbar('Label')      % Initialize/reset and label the bar
%   progressbar(0.5)          % Update
%   progressbar(1)            % Close
%
% Example Function Calls For Multi Bar Usage:
%   progressbar(0, 0)         % Initialize/reset two bars
%   progressbar('A', '')      % Initialize/reset two bars with one label
%   progressbar('', 'B')      % Initialize/reset two bars with one label
%   progressbar('A', 'B')     % Initialize/reset two bars with two labels
%   progressbar(0.3)          % Update 1st bar
%   progressbar(0.3, [])      % Update 1st bar
%   progressbar([], 0.3)      % Update 2nd bar
%   progressbar(0.7, 0.9)     % Update both bars
%   progressbar(1)            % Close
%   progressbar(1, [])        % Close
%   progressbar(1, 0.4)       % Close
%
% Notes:
%   For best results, call progressbar with all zero (or all string) inputs
% before any processing. This sets the proper starting time reference to
% calculate time remaining.
%   Bar color is choosen randomly when the figure is created or reset. Clicking
% the bar will cause a random color change.
%
% Demos:
%     % Single bar
%     m = 500;
%     progressbar % Init single bar
%     for i = 1:m
%       pause(0.01) % Do something important
%       progressbar(i/m) % Update progress bar
%     end
% 
%     % Simple multi bar (update one bar at a time)
%     m = 4;
%     n = 3;
%     p = 100;
%     progressbar(0,0,0) % Init 3 bars
%     for i = 1:m
%         progressbar([],0) % Reset 2nd bar
%         for j = 1:n
%             progressbar([],[],0) % Reset 3rd bar
%             for k = 1:p
%                 pause(0.01) % Do something important
%                 progressbar([],[],k/p) % Update 3rd bar
%             end
%             progressbar([],j/n) % Update 2nd bar
%         end
%         progressbar(i/m) % Update 1st bar
%     end
% 
%     % Fancy multi bar (use labels and update all bars at once)
%     m = 4;
%     n = 3;
%     p = 100;
%     progressbar('Monte Carlo Trials','Simulation','Component') % Init 3 bars
%     for i = 1:m
%         for j = 1:n
%             for k = 1:p
%                 pause(0.01) % Do something important
%                 % Update all bars
%                 frac3 = k/p;
%                 frac2 = ((j-1) + frac3) / n;
%                 frac1 = ((i-1) + frac2) / m;
%                 progressbar(frac1, frac2, frac3)
%             end
%         end
%     end
%
% Author:
%   Steve Hoelzer
%
% Revisions:
% 2002-Feb-27   Created function
% 2002-Mar-19   Updated title text order
% 2002-Apr-11   Use floor instead of round for percentdone
% 2002-Jun-06   Updated for speed using patch (Thanks to waitbar.m)
% 2002-Jun-19   Choose random patch color when a new figure is created
% 2002-Jun-24   Click on bar or axes to choose new random color
% 2002-Jun-27   Calc time left, reset progress bar when fractiondone == 0
% 2002-Jun-28   Remove extraText var, add position var
% 2002-Jul-18   fractiondone input is optional
% 2002-Jul-19   Allow position to specify screen coordinates
% 2002-Jul-22   Clear vars used in color change callback routine
% 2002-Jul-29   Position input is always specified in pixels
% 2002-Sep-09   Change order of title bar text
% 2003-Jun-13   Change 'min' to 'm' because of built in function 'min'
% 2003-Sep-08   Use callback for changing color instead of string
% 2003-Sep-10   Use persistent vars for speed, modify titlebarstr
% 2003-Sep-25   Correct titlebarstr for 0% case
% 2003-Nov-25   Clear all persistent vars when percentdone = 100
% 2004-Jan-22   Cleaner reset process, don't create figure if percentdone = 100
% 2004-Jan-27   Handle incorrect position input
% 2004-Feb-16   Minimum time interval between updates
% 2004-Apr-01   Cleaner process of enforcing minimum time interval
% 2004-Oct-08   Seperate function for timeleftstr, expand to include days
% 2004-Oct-20   Efficient if-else structure for sec2timestr
% 2006-Sep-11   Width is a multiple of height (don't stretch on widescreens)
% 2010-Sep-21   Major overhaul to support multiple bars and add labels
%

persistent progfig progdata lastupdate

% Get inputs
if nargin > 0
    input = varargin;
    ninput = nargin;
else
    % If no inputs, init with a single bar
    input = {0};
    ninput = 1;
end

% If task completed, close figure and clear vars, then exit
if input{1} == 1
    if ishandle(progfig)
        delete(progfig) % Close progress bar
    end
    clear progfig progdata lastupdate % Clear persistent vars
    drawnow
    return
end

% Init reset flag 
resetflag = false;

% Set reset flag if first input is a string
if ischar(input{1})
    resetflag = true;
end

% Set reset flag if all inputs are zero
if input{1} == 0
    % If the quick check above passes, need to check all inputs
    if all([input{:}] == 0) && (length([input{:}]) == ninput)
        resetflag = true;
    end
end

% Set reset flag if more inputs than bars
if ninput > length(progdata)
    resetflag = true;
end

% If reset needed, close figure and forget old data
if resetflag
    if ishandle(progfig)
        delete(progfig) % Close progress bar
    end
    progfig = [];
    progdata = []; % Forget obsolete data
end

% Create new progress bar if needed
if ishandle(progfig)
else % This strange if-else works when progfig is empty (~ishandle() does not)
    
    % Define figure size and axes padding for the single bar case
    height = 0.03;
    width = height * 8;
    hpad = 0.02;
    vpad = 0.25;
    
    % Figure out how many bars to draw
    nbars = max(ninput, length(progdata));
    
    % Adjust figure size and axes padding for number of bars
    heightfactor = (1 - vpad) * nbars + vpad;
    height = height * heightfactor;
    vpad = vpad / heightfactor;
    
    % Initialize progress bar figure
    left = (1 - width) / 2;
    bottom = (1 - height) / 2;
    progfig = figure(...
        'Units', 'normalized',...
        'Position', [left bottom width height],...
        'NumberTitle', 'off',...
        'Resize', 'off',...
        'MenuBar', 'none' );
    
    % Initialize axes, patch, and text for each bar
    left = hpad;
    width = 1 - 2*hpad;
    vpadtotal = vpad * (nbars + 1);
    height = (1 - vpadtotal) / nbars;
    for ndx = 1:nbars
        % Create axes, patch, and text
        bottom = vpad + (vpad + height) * (nbars - ndx);
        progdata(ndx).progaxes = axes( ...
            'Position', [left bottom width height], ...
            'XLim', [0 1], ...
            'YLim', [0 1], ...
            'Box', 'on', ...
            'ytick', [], ...
            'xtick', [] );
        progdata(ndx).progpatch = patch( ...
            'XData', [0 0 0 0], ...
            'YData', [0 0 1 1] );
        progdata(ndx).progtext = text(0.99, 0.5, '', ...
            'HorizontalAlignment', 'Right', ...
            'FontUnits', 'Normalized', ...
            'FontSize', 0.7 );
        progdata(ndx).proglabel = text(0.01, 0.5, '', ...
            'HorizontalAlignment', 'Left', ...
            'FontUnits', 'Normalized', ...
            'FontSize', 0.7 );
        if ischar(input{ndx})
            set(progdata(ndx).proglabel, 'String', input{ndx})
            input{ndx} = 0;
        end
        
        % Set callbacks to change color on mouse click
        set(progdata(ndx).progaxes, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        set(progdata(ndx).progpatch, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        set(progdata(ndx).progtext, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        set(progdata(ndx).proglabel, 'ButtonDownFcn', {@changecolor, progdata(ndx).progpatch})
        
        % Pick a random color for this patch
        changecolor([], [], progdata(ndx).progpatch)
        
        % Set starting time reference
        if ~isfield(progdata(ndx), 'starttime') || isempty(progdata(ndx).starttime)
            progdata(ndx).starttime = clock;
        end
    end
    
    % Set time of last update to ensure a redraw
    lastupdate = clock - 1;
    
end

% Process inputs and update state of progdata
for ndx = 1:ninput
    if ~isempty(input{ndx})
        progdata(ndx).fractiondone = input{ndx};
        progdata(ndx).clock = clock;
    end
end

% Enforce a minimum time interval between graphics updates
myclock = clock;
if abs(myclock(6) - lastupdate(6)) < 0.01 % Could use etime() but this is faster
    return
end

% Update progress patch
for ndx = 1:length(progdata)
    set(progdata(ndx).progpatch, 'XData', ...
        [0, progdata(ndx).fractiondone, progdata(ndx).fractiondone, 0])
end

% Update progress text if there is more than one bar
if length(progdata) > 1
    for ndx = 1:length(progdata)
        set(progdata(ndx).progtext, 'String', ...
            sprintf('%1d%%', floor(100*progdata(ndx).fractiondone)))
    end
end

% Update progress figure title bar
if progdata(1).fractiondone > 0
    runtime = etime(progdata(1).clock, progdata(1).starttime);
    timeleft = runtime / progdata(1).fractiondone - runtime;
    timeleftstr = sec2timestr(timeleft);
    titlebarstr = sprintf('%2d%%    %s remaining', ...
        floor(100*progdata(1).fractiondone), timeleftstr);
else
    titlebarstr = ' 0%';
end
set(progfig, 'Name', titlebarstr)

% Force redraw to show changes
drawnow

% Record time of this update
lastupdate = clock;

end
% ------------------------------------------------------------------------------
function changecolor(h, e, progpatch) %#ok<INUSL>
% Change the color of the progress bar patch

% Prevent color from being too dark or too light
colormin = 1.5;
colormax = 2.8;

thiscolor = rand(1, 3);
while (sum(thiscolor) < colormin) || (sum(thiscolor) > colormax)
    thiscolor = rand(1, 3);
end

set(progpatch, 'FaceColor', thiscolor)

end
% ------------------------------------------------------------------------------
function timestr = sec2timestr(sec)
% Convert a time measurement from seconds into a human readable string.

% Convert seconds to other units
w = floor(sec/604800); % Weeks
sec = sec - w*604800;
d = floor(sec/86400); % Days
sec = sec - d*86400;
h = floor(sec/3600); % Hours
sec = sec - h*3600;
m = floor(sec/60); % Minutes
sec = sec - m*60;
s = floor(sec); % Seconds

% Create time string
if w > 0
    if w > 9
        timestr = sprintf('%d week', w);
    else
        timestr = sprintf('%d week, %d day', w, d);
    end
elseif d > 0
    if d > 9
        timestr = sprintf('%d day', d);
    else
        timestr = sprintf('%d day, %d hr', d, h);
    end
elseif h > 0
    if h > 9
        timestr = sprintf('%d hr', h);
    else
        timestr = sprintf('%d hr, %d min', h, m);
    end
elseif m > 0
    if m > 9
        timestr = sprintf('%d min', m);
    else
        timestr = sprintf('%d min, %d sec', m, s);
    end
else
    timestr = sprintf('%d sec', s);
end
end
