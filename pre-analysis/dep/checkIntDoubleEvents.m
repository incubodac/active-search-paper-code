clases = []
for i = 1:length(EEG.event)
    clases = [clases ; {class(EEG.event(i).duration)}];
end
unique(clases)
losdou= EEG.event(ismember(clases,'double'))
losint= EEG.event(ismember(clases,'int64'))

for i = 1:length(EEG.event)
    if strcmp(EEG.event(i).duration,'int64')
        EEG.event(i).duration = double(EEG.event(i).duration)
    end
end