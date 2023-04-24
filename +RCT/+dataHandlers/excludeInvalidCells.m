function [data] = excludeInvalidCells(data)

validCells = RCT.dataHandlers.getValidCells(data);
locations = [data.NeuronInfo(:).location];
NL = numel(locations);

for ll = 1:NL
    loc = locations(ll);
    data.NeuronInfo(ll).sortingClassification = data.NeuronInfo(ll).sortingClassification(validCells.(loc));
    data.NeuronInfo(ll).neuronID = data.NeuronInfo(ll).neuronID(validCells.(loc));

    for tt = 1:numel(data.trials)
        data.trials(tt).Y.(loc) = data.trials(tt).Y.(loc)(:, validCells.(loc));
    end
end

