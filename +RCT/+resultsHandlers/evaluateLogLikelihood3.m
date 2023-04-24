function [LL_s, LL2_s, FR_s, dim_N_ranges, Z, LL, LL2] = evaluateLogLikelihood3(gmlm, samples, modelSetup, filterLength)
%%
dim_M_train = numel(modelSetup.trialsUsed);
lambda_log_scale = modelSetup.includeTrialwiseLatent_log_scaleFunc(dim_M_train);
log_constant_scale = log(modelSetup.constant_scale);

if(nargin < 4 || isempty(filterLength))
    filterLength = 5;
end

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

%% for each samples
NS = size(samples.W,2);

NS_part = 200;
LL = nan(TT, ceil(NS/NS_part));
FR = nan(TT, ceil(NS/NS_part));
if(filterLength > 1)
    LL2 = nan(TT, ceil(NS/NS_part));
end
LL_part = nan(TT, NS_part);
FR_part = nan(TT, NS_part);
Z = nan(dim_M, NS);

col_counter = 0;
sample_ctr = 0;

for sample_idx = 1:NS
    sample_idx_c = mod(sample_idx - 1, NS_part) + 1;
    if(sample_idx == 1 || sample_idx == NS || mod(sample_idx,50) == 0)
        fprintf("Sample %d / %d (%d)\n", sample_idx, NS, sample_idx_c);
    end


    %% gets the base log rate for the current parameters
    params.W(:) = samples.W(:,sample_idx);
    params.B(:,:) = samples.B(:,:,sample_idx);
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
        LL_part(:,sample_idx_c) = sum(LL_info.logLikeFun(log_rate, Y, dt),2);
        sample_ctr = sample_ctr + 1;

        FR_part(:, sample_idx_c) = mean(exp(log_rate),2);
    else
    
        [log_rate, xx] = gmlm.computeLogRate_host_v2(params);
        T = size(log_rate,1);
        X_stim_1 = xx(1).c(1:T, :, 1);
        X_stim_2 = xx(1).c((1:T) + T, :, 1);
        X_stim_3 = xx(1).c((1:T) + 2*T, :, 1);
        X_stim_4 = xx(1).c((1:T) + 3*T, :, 1);
        X_resp   = xx(2).c(:,:,1);
        X = [X_stim_1, X_stim_2, X_stim_3, X_stim_4, X_resp];
        V = [repmat(params.Groups(1).V, 1, 4) params.Groups(2).V];

        %% gets std of latents
        R_stim = size(params.Groups(1).V,2);
        idx_log_c2 = 1;
        log_c2 = samples.Groups(1).H(idx_log_c2, sample_idx);

        idx_tau = 2;
        log_tau = samples.Groups(1).H(idx_tau, sample_idx);

        idx_phi = 2 + (1:(R_stim));
        log_phi = reshape(samples.Groups(1).H(idx_phi, sample_idx), 1, R_stim);

        idx_lambda = 2 + R_stim + (R_stim*gmlm.dim_P) + (1:(R_stim*12));
        log_lambda_ss = reshape(samples.Groups(1).H(idx_lambda, sample_idx), [], R_stim);
        log_lambda_ss = log_lambda_ss(9:12,:);

        log_unregularized_scale2 = 2*(log_constant_scale + log_phi + log_tau + log_lambda_ss + lambda_log_scale);
        log_scales_stim          = log_unregularized_scale2 + log_c2 - kgmlm.utils.logSumExp_pair(log_c2, log_unregularized_scale2);

        R_resp = size(params.Groups(2).V,2);
        idx_log_c2 = 1;
        log_c2 = samples.Groups(2).H(idx_log_c2, sample_idx);

        idx_tau = 2;
        log_tau = samples.Groups(2).H(idx_tau, sample_idx);

        idx_phi = 2 + (1:(R_resp));
        log_phi = reshape(samples.Groups(2).H(idx_phi, sample_idx), 1, R_resp);

        idx_lambda = 2 + R_resp + (R_resp*gmlm.dim_P) + (1:(R_resp*6));
        log_lambda_ss = reshape(samples.Groups(1).H(idx_lambda, sample_idx), [], R_resp);
        log_lambda_ss = log_lambda_ss(6,:);

        log_unregularized_scale2 = 2*(log_constant_scale + log_phi + log_tau + log_lambda_ss + lambda_log_scale);
        log_scales_resp          = log_unregularized_scale2 + log_c2 - kgmlm.utils.logSumExp_pair(log_c2, log_unregularized_scale2);

        log_scales_c = [reshape(log_scales_stim',[],1); log_scales_resp];
        scales_c = exp(log_scales_c);

        R = 100;

        llc = zeros(TT, R);

        latents_c = randn(numel(scales_c),R) .* scales_c(:); 
        for rr = 1:R
            % REDO THIS SO MONTE CARLO ERRORS AREN'T CORRELATED OVER TRIALS

            llc(:,rr) = sum(LL_info.logLikeFun(X*diag(latents_c(:,rr))*V' + log_rate, Y, dt),2);
        end
        LL_part(:,sample_idx_c) = kgmlm.utils.logMeanExp(llc,2);
        sample_ctr = sample_ctr + 1;


        FR_part(:, sample_idx_c) = mean(exp(log_rate),2);

    end

    %%
    for mm = 1:dim_M
        xx = dim_N_ranges(mm):(dim_N_ranges(mm+1)-1);
        Z(mm, sample_idx) = sum(LL_part(xx,sample_idx_c),1);
    end

    %%
    if(mod(sample_idx_c, NS_part) == 0 || sample_idx == NS)
        col_counter = col_counter + 1;
        FR(:, col_counter) = mean(FR_part(:, 1:sample_idx_c), 2)*(sample_idx_c/NS);
        LL(:, col_counter) = kgmlm.utils.logSumExp(LL_part(:, 1:sample_idx_c), 2) - log(NS);
        if(filterLength > 1)
            LL_c = LL_part(:, 1:sample_idx_c);
            for mm = 1:dim_M
                xx = dim_N_ranges(mm):(dim_N_ranges(mm+1)-1);
                LL_c(xx,:) = conv2(LL_c(xx,:), ones(filterLength,1), "same");
            end
            LL2(:, col_counter) = kgmlm.utils.logSumExp(LL_c, 2) - log(NS);
        end
        %fprintf(" here %d , %d , %d, %d\n", sample_idx, sample_idx_c, col_counter, sample_ctr);
        sample_ctr = 0;
    end
end

FR_s = sum(FR,2);
LL_s = kgmlm.utils.logSumExp(LL(:,1:col_counter),2)  - 0*log(NS);
if(filterLength > 1)
    LL2_s = kgmlm.utils.logSumExp(LL2(:,1:col_counter),2)  - 0*log(NS);
else
    LL2_s = LL_s;
end
Z = kgmlm.utils.logMeanExp(Z, 2);

end
