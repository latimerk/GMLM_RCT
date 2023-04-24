function [validCells] = getValidCells(data)
VALID_CELL_CLASSIFICATIONS = ["good" "mua"];
% VALID_CELL_CLASSIFICATIONS = ["good"];
MIN_FIRING_RATE_SP_S = 5; %spikes per second
WINDOW_MS = 500; % millieseconds (bins)

MIN_SPK_COUNT = MIN_FIRING_RATE_SP_S * (WINDOW_MS / 1e3);



folders = RCT.dataHandlers.checkSession([], data.subject, data.session, [data.NeuronInfo(:).location]);
[~,fname_curation] = RCT.dataHandlers.getProcessedDataFileName(folders, data.subject, data.session);
if(~isempty(fname_curation) && exist(fname_curation, "file"))
    load(fname_curation, "badCells");
    for ff = 1:numel(data.NeuronInfo)
        if(isfield(badCells,data.NeuronInfo(ff).location) && ~isempty( badCells.(data.NeuronInfo(ff).location)))
            bc = ismember(data.NeuronInfo(ff).neuronID, badCells.(data.NeuronInfo(ff).location));
            data.NeuronInfo(ff).sortingClassification(bc) = "relabled_bad";
        end
    end
end

locations = [data.NeuronInfo(:).location];
NL = numel(locations);
NT = numel(data.trials);

validCells = struct();

for ll = 1:NL
    vv = ismember(strrep(data.NeuronInfo(ll).sortingClassification, " ", ""), VALID_CELL_CLASSIFICATIONS);
    vv = vv(:);

    scs = nan(numel(vv), NT);
    for tt = 1:NT
        ww = data.trials(tt).noise_on + (0:(WINDOW_MS-1));
        scs(:, tt) = sum(data.trials(tt).Y.(locations(ll))(ww, :), 1)';
    end
    scs = mean(scs,2);

    vv = vv & scs >= MIN_SPK_COUNT;
    validCells.(locations(ll)) = vv;
end


