function [data, session_name] = loadData(subject, session, locations, removeInvalidCells, downsample)
if(nargin < 3)
    locations = [];
end
if(nargin < 4 || isempty(removeInvalidCells))
    removeInvalidCells = true;
end
if(nargin < 5 || isempty(downsample))
    downsample = -1;
end

[folders, subject, ~, session_name, locations] = RCT.dataHandlers.checkSession([], subject, session, locations);

fname = RCT.dataHandlers.getProcessedDataFileName(folders, subject, session_name);
if(~exist(fname, "File"))
    fprintf("Processed data not found. Processing... ");
    data = RCT.dataHandlers.processRawData(subject, session_name);
    fprintf("done.\n");
else
    data = load(fname, "data");
    data = data.data;
end

all_locations = [data.NeuronInfo(:).location];

vv = ismember(all_locations, locations);
data.NeuronInfo = data.NeuronInfo(vv);
for tt = 1:numel(data.trials)
    fs = cellfun(@string,fieldnames(data.trials(tt).Y));
    for ff = 1:numel(fs)
        if(~ismember(fs(ff), locations))
            data.trials(tt).Y = rmfield(data.trials(tt).Y, fs(ff));
        end
    end
end



if(removeInvalidCells)
    data = RCT.dataHandlers.excludeInvalidCells(data);
end

%downsample
if(downsample > 1)
    downsample = ceil(downsample);
    data.bin_size_ms = downsample * data.bin_size_ms;

    fields_to_fix = ["fix_on", "noise_on", "targets_on", "fix_off", "stim_on", "saccade_end"];
    for tt = 1:numel(data.trials)
        fs = cellfun(@string,fieldnames(data.trials(tt).Y));
        for ff = 1:numel(fs)
            Y = data.trials(tt).Y.(fs(ff));
            Y = Y(1:(floor(size(Y,1)/downsample)*downsample), :);

            data.trials(tt).Y.(fs(ff)) = squeeze(sum(reshape(Y,[downsample size(Y,1)/downsample size(Y,2)]),1));
        end

        for ff = 1:numel(fields_to_fix)
            data.trials(tt).(fields_to_fix(ff)) = ceil(data.trials(tt).(fields_to_fix(ff)) / downsample);
        end
    end

end