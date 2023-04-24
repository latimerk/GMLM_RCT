function [data_raw, session_name] = loadRawData(folders, subject, session, locations)
if(nargin < 4)
    locations = [];
end

[folders, subject, ~, session_name, locations] = RCT.dataHandlers.checkSession(folders, subject, session, locations);

data_raw = struct("Bhv", cell(numel(locations),1), "Neuro", [], "location", []);
for ii = 1:numel(locations)
    fname = RCT.dataHandlers.getRawDataFileName(folders, subject, session_name, locations(ii));

    for jj = 1:numel(fname)
        if(exist(fname(jj), "file"))
            load(fname(jj), "data");
            data_raw(ii).Bhv = data.Bhv;
            data_raw(ii).Neuro = data.Neuro;
            data_raw(ii).location = locations(ii);
            break;
        end
    end
end