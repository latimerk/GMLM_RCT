function [LL, LL_prev, LL_curr, SpikeRate] = evaluateLogLikelihood(gmlm, samples, modelSetup)



%%
dim_M = gmlm.dim_M;
dim_M_train = numel(modelSetup.trialsUsed);
lambda_log_scale = log(modelSetup.includeTrialwiseLatent_scaleFunc(dim_M_train));
log_constant_scale = log(modelSetup.constant_scale);

gmlm.setupComputeStructuresHost();
for jj = 1:numel(samples.Groups)
    gmlm.setDimR(jj, size(samples.Groups(jj).V,2));
end

params = gmlm.getEmptyParamStruct();
LL_info = gmlm.LL_info;
Y = LL_info.Y;
dt = LL_info.bin_size;
dim_N_ranges = LL_info.dim_N_ranges;

opts = optimoptions("fminunc", "algorithm", "quasi-newton", "gradobj", "on", "hessian", "off", "display", "off");

%% for each samples
NS = size(samples.W,2);


SpikeRate = [];%cell(gmlm.dim_M,1);
LL = cell(gmlm.dim_M,1);
LL_prev = [];%cell(gmlm.dim_M,1);
LL_curr = [];%cell(gmlm.dim_M,1);
P = gmlm.dim_P;
for mm = 1:numel(LL)
    nn = gmlm.dim_N(mm);
    %SpikeRate{mm} = zeros(nn, P, NS,2);
    LL{mm} = zeros(nn, P, NS);
    %LL_prev{mm} = zeros(nn, NS);
    %LL_curr{mm} = zeros(nn, NS);
end

ff = exp(-(-200:200).^2./(2*15^2));
ff = ff(:)./sum(ff)*1e3;

for sample_idx = 1:NS
    if(sample_idx == 1 || sample_idx == NS || mod(sample_idx,50) == 0)
        fprintf("Sample %d / %d\n", sample_idx, NS);
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


    if(~modelSetup.includeTrialwiseLatent)
        log_rate = gmlm.computeLogRate_host_v2(params);
        %% 
        ll = LL_info.logLikeFun(log_rate, Y, dt);

        %% for each trial
        for mm = 1:gmlm.dim_M
    
            tts = dim_N_ranges(mm):(dim_N_ranges(mm+1)-1);
            LL{mm}(:,:,sample_idx) = ll(tts,:);
            %SpikeRate{mm}(:,:,sample_idx,1) = exp(log_rate(tts,:));
            %SpikeRate{mm}(:,:,sample_idx,2) = conv2(Y(tts,:),ff,"same");
    
        end
    else
    
        [log_rate, xx] = gmlm.computeLogRate_host_v2(params);
        T = size(log_rate,2);
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

        %% for each trial
        for mm = 1:dim_M
            if(mm == 1 || mm == dim_M || mod(mm,30) == 0)
                fprintf("\tTrial %d / %d\n", mm, dim_M);
            end

            tts = dim_N_ranges(mm):(dim_N_ranges(mm+1)-1);
        
            %% build design matrix
            c_trial = log_rate(tts,:);
            X_trial = X(tts,:);
            Y_trial = Y(tts,:);

            %% predict ahead for each time point: log p(y_{t, p} | theta, x, y_{1:t-1, :})
            % for first time point, only have the prior
            mu_0   = zeros(numel(scales_c),1);
            logdet_0  = sum(log_scales_c);
            invsig2_0 = diag(scales_c); 
            z_p = zeros(size(X_trial));

            invsig2_mu_0 = zeros(numel(scales_c),1);

            for tt = 1:size(X_trial,1)
                %% maximize for latents and do Gaussian approximation
                X_c = X_trial(1:tt,:);
                c_c = c_trial(1:tt,:);
                Y_c = Y_trial(1:tt,:);

                warning("CHECK DERIVATIVES!");
                nlp = @(w)nlpFunc(w, X_c, V, scales_c,  @(rr )LL_info.logLikeFun(rr + c_c, Y_c, dt));
                mu_c = fminunc(nlp, mu_0, opts);

                [~,~,invsig2_c] = nlp(mu_c); 

                logdet_c = kgmlm.utils.logdet(invsig2_c);

                invsig2_mu_c = invsig2_c*mu_c;

                %z_c = (mu_c + mu_0)./2;
                z_c = (invsig2_0 + invsig2_c)\(invsig2_mu_0 + invsig2_mu_c);
                z_p(tt,:) = z_c;

                LL_prev{mm}(tt,sample_idx) = -1/2*((z_c - mu_0)'*(invsig2_0*(z_c - mu_0))) + 1/2*logdet_0;
                LL_curr{mm}(tt,sample_idx) = -1/2*((z_c - mu_c)'*(invsig2_c*(z_c - mu_c))) + 1/2*logdet_c;

                %% set prev time point
                mu_0 = mu_c;
                invsig2_0 = invsig_c;
                logdet_0 = logdet_c;
                invsig2_mu_0 = invsig2_mu_c;
            end
            
            LL{mm}(:,:,sample_idx) = LL_info.logLikeFun((X_trial.*z_p)*V' + c_trial, Y_trial, dt);
            error("NEED TO FIX THIS!")
        end
    end
end

end


function [f,g,h] = nlpFunc(w, X, V, lambda, llFunc)

log_rate = (X.*w(:)')*V';

if(nargin > 2)
    [lls, dlls, d2lls] = llFunc(log_rate); 
elseif(nargin > 1)
    [lls, dlls] = llFunc(log_rate); 
else
    lls = llFunc(log_rate); 
end

f = -sum(lls, "all") + 1/2*sum(w.^2.*lambda,"all");
if(nargin > 1)
    g = w.*lambda - sum((X'*dlls).*V',2);
end
if(nargin > 2)
    h = diag(lambda);
    for pp = 1:size(V,1)
        X2 = (X).*d2lls(:,pp);
        h = h + (X2'*X2).*(V(:,pp)*V(:,pp)');
    end
end

end