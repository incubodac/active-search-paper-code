function [WinRej] = bad_ET_ArtifactDetect(EEG,WinRej)
%after runnin function uf_continuousArtifactDetect() function from 
%unfold toolbox this function adds new windows to exclude from
%adeconvolution analysis. Events to exclude are marked as 'bad_ET' 
% input:
% EEG: EEGLAB struct
% Winrej: nx2 matrix output of uf_continuousArtifactDetect
%from c.r.a.p. detection 
%utput:
%Winrej: same as input with new windows to rejct coming from 'bad_ET'
%events. Overlaped intervals were not corrected so there will be some
%redundant info.
Nsamples = length(EEG.data);
badEvents =  strcmp({EEG.event.type},'bad_ET');
Wintmp =[];

if ~isempty([EEG.event(badEvents).latency])
    Wintmp(:,1) = [EEG.event(badEvents).latency]';
    durtmp      = [EEG.event(badEvents).duration]';
    Wintmp(:,2) = Wintmp(:,1)+ durtmp ;
    if Wintmp(end,2) > Nsamples
        Wintmp(end,2) = Nsamples;
    end
end

WinRej = [WinRej  ; Wintmp];


