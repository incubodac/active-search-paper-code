classdef FP_dac
    properties
        cfg 
    end
    methods
        function obj = FP_dac()
            fprintf('----------------------------------------------\n')
            fprintf('\t\trepo: FP_dac: prehybridana\n')
            fprintf('----------------------------------------------\n')            
            %general use for passinginput and ouput files to pre-processing fucntions,
            obj.cfg.inFile               = [];
            obj.cfg.outfile              = [];
            %parameters
            obj.cfg.eye                  = 'R';%eye to be used 'R' or 'L'
            obj.cfg.ref                  = [129 130];%vector of references for pop_biosig
            obj.cfg.nChans               = 128; %number of chans with data
            obj.cfg.screenSize           = [1920 1080]; %monitor resolution [x y]
            %obj.cfg.chanlocsFilePath     = '/functions/resources/Standard-10-5-Cap385.sfp'; %relative chanloc file               
            %ETEEG.event
            obj.cfg.keyword              = 'b\''ETSYNC';%keyword to parse eyetracker data
            %EEG prefiltering
            obj.cfg.lineFreq             =  50; 
            obj.cfg.preFilterEdges       = [0.1 100];%edges for band pass filter using pop_eegfiltnew()

            obj.cfg.noFiltAnalog         = 0; %[Bool] wheter to remove analog before filtering (1) or not (0)
            %interpolation
            obj.cfg.nHeadChans           = 128; %Number of channels to inspect  before interpolation
            obj.cfg.badchanslist         = []; %cell with {'subj', [badchans]} update  after inspetion
            % Load EEG and ET files and syn using pop_importeyetracker
            obj.cfg.inFileEEG            = []; %full path to .set EEG file
            obj.cfg.inFileET             = []; %full path to .mat ET file
            obj.cfg.marks                = [150 200];%marks to sync ET with EEG   
            % ICA training
            obj.cfg.icaFilterEdge        = 2; %Hz
            obj.cfg.trialDur             = 4; %seg
            obj.cfg.preDur               = 6; %seg 
            obj.cfg.stimMark             = '255';
            obj.cfg.preMark              = '151';
            %obj.cfg.preTrialMark         = '12';
            %obj.cfg.posTrialMark         = '14';
            % Epoched data
            obj.cfg.postFilterEdges      = [0.2 40];
            % IC's selection
            obj.cfg.weights             = []; %full path to .mat file with weights struct
            obj.cfg.badcomp             = []; %channels to remove based on Plöchl
        end
    end
    methods(Static)
    function names = loadSubjects(folder, expe, pre, type, keyw)     %Names of files to load (pre and expe are only used when data were recorded separately)
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
    function ET  = renameEtEvents(fileIn, fileOut, keyword)    %asc file path as in and .mat struct out to parse ET data
        fprintf('DAC (2022-05-06): renameEtEvents: remove all events except 200 255\n')
        ET = parseeyelink(fileIn, fileOut, keyword);  
        %keyboard
        %parsing ending event
%         test_end         = regexp(ET.messages,['MSG\s(\d+)\s+','KEY','\s?(\s+)'] ,'tokens')';
%         test_end         = [test_end{:}]';
%         test_end         = cellfun(@str2double,test_end,'UniformOutput',false);

        %parsing eyemap events
%         test             = regexp(ET.messages,['MSG\s(\d+) b\''ETSYNC ' keyword '\'''],'tokens')';
%         test2            = [test{:}]';
%         test3            = cellfun(@str2double,test2,'UniformOutput',false);
%         test3 = cellfun(@(x) [test3(1) 151] ,test3,'UniformOutput',0);
%         test3            = cellfun(@str2double,test3,'UniformOutput',false);
% 
%         %keyboard
%         %6/02/20 try using test_end events for eye map detetion too 
% %         end_evt          = cleantriggers(cell2mat(test_end));       %trigger and rising flank of trigger for END's events
%         eyemap_evt       = cleantriggers(cell2mat(test3));%CHANGE HERE OLD WITH     test3          %trigger and rising flank of trigger for eyemap
% 
% 
%         ET.event         =  eyemap_evt;
% 
%         idxcelleyeall_2  = ismember(ET.event(:,1),ET.event(:,1));%?
%         figure;plot(diff(cell2mat({ET.event(idxcelleyeall_2,1)})/1000),'b.--')% time difference between events ET
% 
% 
%         tol=0.05;
%         ind_EMAP_long    = find(abs(diff(cell2mat({ET.event(idxcelleyeall_2,1)})/1000)-6)<tol);
%         figure;plot(ind_EMAP_long,'b.--');title('Check Eyemap Parsing');
%         fprintf('Modifying events in ET')
%         event_times_EMAP_eye  =  ET.event(ind_EMAP_long,1);
% 
% 
%         %EYEMAP EYE events: long 6sec gap 
%         ET.event_old     =  ET.event; %save back-up of original events
%         ET.event         =  ET.event(ind_EMAP_long)
% 
%         for ii=1:numel(event_times_EMAP_eye)
% 
%             ET.event(ii,1)   = event_times_EMAP_eye(ii);
%             ET.event(ii,2)   = 290;
%         end
%         ET.event(ii+1,1)    = end_evt(end,1); %add last event 
%         ET.event(ii+1,2)    = 300; 
%         fprintf('Overwriting ET struct, new events added')
%         %output file must be complete path and .mat file 
          ET.event(find(ET.event(:,2) == 50),:)     =  [];
          ET.event(find(ET.event(:,2) == 150),:)    =  [];
          ET.event(find(ET.event(:,2) == 151),:)    =  [];
          ET.event(find(ET.event(:,2) == 200),:)    =  [];
          ET.event(find(ET.event(:,2) == 201),:)    =  [];
          ET.event(find(ET.event(:,2) == 251),:)    =  [];
          
          ET.event(find(ET.event(:,2) == 250),2)    = 255;
          ET.event(find(ET.event(:,2) == 255,1,'last'),2) = 200;
          save(fileOut,'-struct','ET')% overwrites mat file the same way parseeyelink does it.
    end
    function EEG = renameEEGevents(EEG)
                tol                =  0.05;
                idxcellall         =  true(1,numel(EEG.event));
                ind_EMAP_longEEG   =  find( abs(diff(cell2mat({EEG.event(idxcellall).latency})/EEG.srate)-6)<tol ); % Asi esta tomando el primero de los dos de la diferencia. Si no hay que sumarle 1
                events_EMAP_EEG    =  cell2mat({EEG.event(ind_EMAP_longEEG).latency})';
                %EYEMAP EEG events: long 6 sec gap
                %MOMENTARILY I'M NOT CHECKING BEFORE SYNCHRONIZE
                %figure;plot(event_times_EMAP_eye,events_EMAP_EEG,'go--')
                %[rho,pval]         =corr(event_times_EMAP_eye,events_EMAP_EEG)


                fprintf('Modifico eventos en EEG')
%                 %defines marks for eyemap events    
%                 for ii=1:numel(ind_EMAP_longEEG)
%                     EEG.event( ind_EMAP_longEEG(ii) ).type = 290;
%                 end
%                 %find last event of experiment
%                 indlast = find(ismember([EEG.event.type],[180]));
% 
%                 EEG.event(indlast(end)).type = 300
%                 %here EEG struct should have the correct marks for
                %synchronization.
                EEG              = eeg_checkset( EEG );
    end
    function EEG = addChanInfo(EEG)
      if EEG.nbchan == 71
        EEG = pop_select( EEG,'nochannel',65:71);
        EEG = eeg_checkset( EEG );
       elseif EEG.nbchan == 72
 
       EEG = pop_select( EEG,'nochannel',65:72);
        EEG = eeg_checkset( EEG );
      end
    end

    function EEG = addFakeEOG(EEG)
        
        
    
    end 
    function [EEG badcahns] = prepareDataForICA(inFile, icaFilterEdge,cfg)      
    % DAMIAN
          % pasar todo esto  a tu funcion custom para epochear                                      
        trialDur         = cfg.trialDur;
        preDur           = cfg.preDur;
        stimMark         = cfg.stimMark; 
        preMark          = cfg.preMark; 
        EEG                          = pop_loadset(inFile);
        EEG                          = eeg_checkset(EEG);  
        EEG.info.filter.loicatrain   = icaFilterEdge;
        
        EEG =   pop_eegfiltnew(EEG, EEG.info.filter.loicatrain,[]);      
        

        if cfg.epoch_before_ica
            % Epoch del pre - marcas de pre [-0.2 preDur]
            %keyboard
            EEGpre = pop_epoch(EEG,{preMark},[-0.2 preDur]);
    
            
            % Epoch de estimulos - Inicio del trial [-0.2 TrialDur]
          %  EEGstim = pop_epoch(EEG,{stimMark},[-0.2 trialDur]); DAC 2022 may
          %  17
    
            
            % merge pre y estimulos
            %EEG = pop_mergeset(EEGpre,EEGstim) this is not possible bc
            %differences in length
            
            EEGpre  = eeg_epoch2continuous(EEGpre)
            %EEGstim = eeg_epoch2continuous(EEGstim) DAC may 17
            %%%%%%%%%%%%%%%%%%%%%DAC may 17
            Eind = find(strcmp({EEG.event.type},'151'));
            Sind = find(strcmp({EEG.event.type},'255'));
    
            eyemaps_start = [EEG.event([1:5:31]).latency]; %%DAC bug detected may24 correct after may
            stims_start = [EEG.event([1:30:210]).latency]; %%DAC bug detected may24 correct after may
    
            eyemaps_ranges = [eyemaps_start(1) stims_start(1);...
                              eyemaps_start(2) stims_start(2);...
                              eyemaps_start(3) stims_start(3);...
                              eyemaps_start(4) stims_start(4);...
                              eyemaps_start(5) stims_start(5);...
                              eyemaps_start(6) stims_start(6);...
                              eyemaps_start(7) stims_start(7)];
            %%%%%%%%%%%%%%%%%%%%%%
            %keyboard
            EEGstim = pop_select(EEG, 'nopoint', eyemaps_ranges,'channel', 1:148)
            %%%%%%%%%%%%%%%%%%%%%%%
            EEG = pop_mergeset(EEGpre,EEGstim)      
        end
    end
    end%methods
end% classdef
