function [params] = generateInitialParameters(gmlm, modelSetup)

params = gmlm.getRandomParamStruct();
if(~isfield(modelSetup, "useRidge") || ~modelSetup.useRidge)
    params.B(:) = randn(gmlm.dim_B, gmlm.dim_P)*1;
    params.H(:) = randn(size(params.H))*0.1;
    params.H(2) = randn*0.1;
    params.H((2 + gmlm.dim_B + 1):end) = randn(numel(params.H) - 2 - gmlm.dim_B,1)*0.1;

    A = mean(cell2mat({gmlm.trials(:).Y}'));
    
    mu = log(A) - log(gmlm.bin_size);


%     params.H(1) = mean(params.W);
%     params.W(:) = randn(gmlm.dim_P,1)*1;
    params.H(1) = mean(mu) + randn*0.1;
    params.W(:) = (mu(:) - mean(mu) + randn(gmlm.dim_P,1)*0.1)./exp(params.H(2));
end

for jj = 1:numel(params.Groups)
    if(strcmpi(params.Groups(jj).name, "Trial_Variability"))
        params.Groups(jj).T{end}(:,:) = randn(size(params.Groups(jj).T{end}))*1e-2;
    end
end

includeFixation = isfield(modelSetup, "fixationConfig") && ~isempty(modelSetup.fixationConfig);
for jj = (3+includeFixation):numel(params.Groups)
    params.Groups(jj).T{2}(:,:) = randn(size(params.Groups(jj).T{2})) * 1e-2;
    for ss = 3:numel(params.Groups(jj).T)
        params.Groups(jj).T{ss}(1,:) = 1;
    end
end