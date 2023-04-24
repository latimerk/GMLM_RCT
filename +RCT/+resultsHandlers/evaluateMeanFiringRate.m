function [FR_s, dim_N_ranges] = evaluateMeanFiringRate(gmlm, samples, modelSetup) %, lambda_s
%%
%dim_M_train = numel(modelSetup.trialsUsed);
%lambda_log_scale = log(modelSetup.includeTrialwiseLatent_scaleFunc(dim_M_train));
%log_constant_scale = log(modelSetup.constant_scale);

gmlm.setupComputeStructuresHost();
for jj = 1:numel(samples.Groups)
    gmlm.setDimR(jj, size(samples.Groups(jj).V,2));
end

params = gmlm.getEmptyParamStruct();
LL_info = gmlm.LL_info;
Y = LL_info.Y;
dt = LL_info.bin_size;
dim_N_ranges = LL_info.dim_N_ranges;
dim_M = numel(dim_N_ranges) - 1;
TT = size(Y,1);

dim_B = gmlm.dim_B();

%% for each samples
NS = size(samples.W,2);

P = gmlm.dim_P;

NS_part = 100;
FR = nan(TT, P, ceil(NS/NS_part));
FR_part = nan(TT, P, NS_part);
% lambda      = nan(TT, P, ceil(NS/NS_part));
% lambda_part = nan(TT, P, NS_part);

col_counter = 0;
sample_ctr = 0;

for sample_idx = 1:NS
    sample_idx_c = mod(sample_idx - 1, NS_part) + 1;
    if(sample_idx == 1 || sample_idx == NS || mod(sample_idx,50) == 0)
        fprintf("Sample %d / %d (%d)\n", sample_idx, NS, sample_idx_c);
    end

    %% gets the base log rate for the current parameters
    params.W(:) = samples.W(:,sample_idx);
    if(dim_B > 0)
        params.B(:,:) = samples.B(:,:,sample_idx);
    end
    for jj = 1:numel(samples.Groups)
        params.Groups(jj).V(:,:) = samples.Groups(jj).V(:,:,sample_idx);
        for ss = 1:numel(params.Groups(jj).T)
            params.Groups(jj).T{ss}(:,:) = samples.Groups(jj).T{ss}(1:size(params.Groups(jj).T{ss},1),:,sample_idx);
        end
    end

    %%
    if(~modelSetup.includeTrialwiseLatent)
        log_rate = gmlm.computeLogRate_host_v2(params);
        %% 
        FR_part(:, :, sample_idx_c) = 1 - exp(-exp(log_rate + log(dt))); % probability of spike under the truncated Poisson model
%         lambda_part(:, :, sample_idx_c) = log_rate; 
        sample_ctr = sample_ctr + 1;
    else
        error("Not implemented");
    end

    %%
    if(mod(sample_idx_c, NS_part) == 0 || sample_idx == NS)
        col_counter = col_counter + 1;
        FR(:, :, col_counter) = mean(FR_part(:, :, 1:sample_idx_c), 3)*(sample_idx_c/NS);
%         lambda(:, :, col_counter) = mean(lambda_part(:, :, 1:sample_idx_c), 3)*(sample_idx_c/NS);
        %fprintf(" here %d , %d , %d, %d\n", sample_idx, sample_idx_c, col_counter, sample_ctr);
        sample_ctr = 0;
    end
end

FR_s = sum(FR,3);
% lambda_s = sum(lambda,3);

end
