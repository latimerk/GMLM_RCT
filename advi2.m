NI = 6000;
% alphaSteps.step  = [1    100  200   2000   5000];
% alphaSteps.alpha = [1e-4 1e-3 1e-2  1e-3   1e-4];
alphaSteps.step  = [1    3000   5000];
alphaSteps.alpha = [1e-2 1e-3   1e-4];
alphaSteps.reset = 1;%[1 1000 2000 3000];
M = 10;

% ADVI/ BBVI: full rank

opts = gmlm.getComputeOptionsStruct(true);
res = gmlm.getEmptyResultsStruct(opts);
w_init = gmlm.vectorizeParams(params, opts);
P = numel(w_init);
nlpostFunction = @(ww) gmlm.vectorizedNLPost_func(ww, params, opts, res);
nllFunction = @(ww) gmlm.vectorizedNLL_func(ww, params, opts, res);

params_0 = params;
params_0.B(:) = 0;
params_0.W(:) = 0;
params_0.H(:) = 0;
for jj = 1:numel(params_0.Groups)
    params_0.Groups(jj).V(:) = 0;
    params_0.Groups(jj).T{1}(:) = 0;
    params_0.Groups(jj).T{2}(:) = 0;
    params_0.Groups(jj).H(:) = 0;
    params_0.Groups(jj).H(1) = 1;
end
c_mask = gmlm.vectorizeParams(params_0, opts) > 0;


scaled_WB = isfield(gmlm.GMLMstructure, "scaleParams") && ~isempty(gmlm.GMLMstructure.scaleParams);
J = numel(params.Groups);
scaled_VT = false(J,1);
for jj = 1:J
    if(isfield(gmlm.GMLMstructure.Groups(jj), "scaleParams") && ~isempty(gmlm.GMLMstructure.Groups(jj).scaleParams))
        scaled_VT(jj) = true;
    end
end

mu = w_init;
L  = ones(P,1)*log(0.1);

mu(c_mask) = log(10);
L(c_mask)  = log(0.1);
c_mask_iters = 1000;

elp = nan(NI, 1);
ell = nan(NI, 1);
H = nan(NI, 1);
nl_mu = nan(NI, 1);
nl_mu_prev = -nllFunction(mu);
ell_prev = nl_mu_prev;



for ii = 1:NI
    fprintf("ii = %d / %d (%e, %e)\n", ii, NI,nl_mu_prev,ell_prev);
    eta  = randn(P, M);
    zeta = exp(L).*eta + mu;

    ll = nan(M,1);
    lp = nan(M,1);
    dlp = nan(P, M);

    for mm = 1:M
        [nlpost, ndW, ~, results] = nlpostFunction(zeta(:,mm));
        lp(mm) = -nlpost;
        ll(mm) = results.log_likelihood;
        dlp(:,mm) = -ndW;
    end
    params1 = gmlm.devectorizeParams(zeta(:,mm), params, opts);
    params2 = gmlm.rescaleParamStruct(params1, scaled_WB, scaled_VT);
    RCT.utils.printRegularizedParameterStatus(params2);

    mu_step = mean(dlp,2);
    L_step  = mean(dlp.*eta,2).*exp(L) + 1;


%     tau = 1;
%     beta_mu = 0.01; % eta in paper
%     beta_L = 0.1; 
%     alpha = 0.1;
%     epsilon = 1e-16;
%     if(ii == 1)
%         s_mu = mu_step.^2;
%         s_L = L_step.^2;
%     else
%         s_mu = alpha*mu_step.^2 + (1-alpha).*s_mu;
%         s_L = alpha*L_step.^2 + (1-alpha).*s_L;
%     end
%     p_mu = beta_mu*ii^(-1/2+epsilon)*(tau + sqrt(s_mu)).^-1;
%     p_L  = beta_L*ii^(-1/2+epsilon)*(tau + sqrt(s_L)).^-1;
% 
%     mu = mu + p_mu.*mu_step;
%     L = L + p_L.*L_step;

    if(ismember(ii, alphaSteps.reset))
        m_mu = zeros(size(mu));
        v_mu = zeros(size(mu));
        m_L = zeros(size(L));
        v_L = zeros(size(L));
    end


    if(ismember(ii, alphaSteps.step))
        [~,jj] = find(alphaSteps.step == ii,1,"first");
        a_c = alphaSteps.alpha(jj);
        alpha_mu = a_c;
        alpha_L = a_c;
    end
    beta_1 = 0.9;
    beta_2 = 0.99;
    epsilon = 1e-8;
    m_mu = beta_1 * m_mu + (1-beta_1)*mu_step;
    v_mu = beta_2 * v_mu + (1-beta_2)*mu_step.^2;
    m_hat_mu = m_mu./(1-beta_1^ii);
    v_hat_mu = v_mu./(1-beta_2^ii);
    mu = mu + alpha_mu*m_hat_mu./(sqrt(v_hat_mu) + epsilon);

    if(ii > 0)
        m_L = beta_1 * m_L + (1-beta_1)*L_step;
        v_L = beta_2 * v_L + (1-beta_2)*L_step.^2;
        m_hat_L = m_L./(1-beta_1^ii);
        v_hat_L = v_L./(1-beta_2^ii);
        L = L + alpha_L*m_hat_L./(sqrt(v_hat_L) + epsilon);
    end

    if(ii < c_mask_iters)
        mu(c_mask) = log(10);
        L(c_mask)  = log(0.1);
        m_mu(c_mask) = 0;
        v_mu(c_mask) = 0;
        m_L(c_mask) = 0;
        v_L(c_mask) = 0;
    end

    nl_mu_prev = -nllFunction(mu);
    nl_mu(ii) = nl_mu_prev;

    elp(ii) = mean(lp);
    ell(ii) = mean(ll);
    ell_prev = ell(ii);
    %(ii) = sum(log(diag(L)));

    if(mod(ii,50) == 0)
        %%
        kgmlm.utils.sfigure(2);
        clf
        subplot(1,3,1)
        
%         if(ii > 3000) 
%             plot(3000:NI, [nl_mu(3000:end) ell(3000:end)]);
%         elseif(ii > 2000) 
%             plot(2000:NI, [nl_mu(2000:end) ell(2000:end)]);
%         elseif(ii > 1000) 
%             plot(1000:NI, [nl_mu(1000:end) ell(1000:end)]);
%         elseif(ii > 500) 
%             plot(500:NI, [nl_mu(500:end) ell(500:end)]);
%         elseif(ii > 200) 
%             plot(200:NI, [nl_mu(200:end) ell(200:end)]);
%         else
%             plot([nl_mu ell]);
%         end
        if(ii > 500) 
            plot(300:NI, [nl_mu(300:end) ell(300:end)]);
        else
            plot([nl_mu ell]);
        end
        subplot(1,3,2)
        plot(mu)
        subplot(1,3,3)
        plot(exp(L))
        drawnow;
        pause(1);
    end
end

params_advi_0 = gmlm.devectorizeParams(mu, params, opts);
params_advi_r = gmlm.rescaleParamStruct(params_advi_0, scaled_WB, scaled_VT);
% clear;
%%
NS = 1000;
samples_advi_r = params_advi_r;
samples_advi_r.W = zeros(size(samples_advi_r.W,1), NS);
samples_advi_r.B = zeros(size(samples_advi_r.B,1), size(samples_advi_r.B,2), NS);
for jj = 1:numel(samples_advi_r.Groups)
    samples_advi_r.Groups(jj).V = zeros(size(samples_advi_r.Groups(jj).V,1), size(samples_advi_r.Groups(jj).V,2), NS);
    samples_advi_r.Groups(jj).T{1} = zeros(size(samples_advi_r.Groups(jj).T{1},1), size(samples_advi_r.Groups(jj).T{1},2), NS);
    samples_advi_r.Groups(jj).T{2} = zeros(size(samples_advi_r.Groups(jj).T{2},1), size(samples_advi_r.Groups(jj).T{2},2), NS);
end

eta  = randn(P, NS);
zeta = exp(L).*eta + mu;

for mm = 1:NS
    params1 = gmlm.devectorizeParams(zeta(:,mm), params, opts);
    params2 = gmlm.rescaleParamStruct(params1, scaled_WB, scaled_VT);

    samples_advi_r.W(:,mm) = params2.W;
    samples_advi_r.B(:,:,mm) = params2.B;
    for jj = 1:numel(samples_advi_r.Groups)
        samples_advi_r.Groups(jj).V(:,:,mm) = params2.Groups(jj).V;
        samples_advi_r.Groups(jj).T{1}(:,:,mm) = params2.Groups(jj).T{1};
        samples_advi_r.Groups(jj).T{2}(:,:,mm) = params2.Groups(jj).T{2};
    end
end

save("adviRun_m.mat", "-v7.3", "params_advi_0", "params_advi_r", "mu", "L", "samples_advi_r");

