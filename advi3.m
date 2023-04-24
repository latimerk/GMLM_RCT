NI = 20000;
M = 100;

% ADVI/ BBVI: full rank

opts = gmlm.getComputeOptionsStruct(true);
res = gmlm.getEmptyResultsStruct(opts);
w_init = gmlm.vectorizeParams(params, opts);
P = numel(w_init);
nlpostFunction = @(ww) gmlm.vectorizedNLPost_func(ww, params, opts, res);
nllFunction = @(ww) gmlm.vectorizedNLL_func(ww, params, opts, res);

scaled_WB = isfield(gmlm.GMLMstructure, "scaleParams") && ~isempty(gmlm.GMLMstructure.scaleParams);
J = numel(params.Groups);
scaled_VT = false(J,1);
for jj = 1:J
    if(isfield(gmlm.GMLMstructure.Groups(jj), "scaleParams") && ~isempty(gmlm.GMLMstructure.Groups(jj).scaleParams))
        scaled_VT(jj) = true;
    end
end

% mu_f = w_init;
% L_f  = 0.01*eye(P);
% mu_f = mu;
% L_f = diag(exp(L));

lls = tril(ones(size(L_f)),-1) == 1;
lls = lls(:);

NN = numel(mu_f);
bb = ones(NN, NN);
bb = triu(bb, -512);
bb = bb ~= 1;
bb = tril(bb);
bb = bb(:);


elp = nan(NI, 1);
ell = nan(NI, 1);
H = nan(NI, 1);
nl_mu = nan(NI, 1);
nl_mu_prev = -nllFunction(mu_f);
ell_prev = nl_mu_prev;

% alpha = 0.0001;
% beta_1 = 0.99;
alpha = 0.001;
beta_1 = 0.9;
beta_2 = 0.999;
epsilon = 10e-8;
for ii = 1:NI
    fprintf("ii = %d / %d (%e, %e)\n", ii, NI,nl_mu_prev,ell_prev);
    eta  = randn(P, M);
    zeta = L_f*eta + mu_f;

    ll = nan(M,1);
    lp = nan(M,1);
    dlp = nan(P, M);

    for mm = 1:M
        [nlpost, ndW, ~, results] = nlpostFunction(zeta(:,mm));
        lp(mm) = results.log_post;
        ll(mm) = results.log_likelihood;
        dlp(:,mm) = -ndW;
    end
    params1 = gmlm.devectorizeParams(zeta(:,mm), params, opts);
    params2 = gmlm.rescaleParamStruct(params1, scaled_WB, scaled_VT);
    RCT.utils.printRegularizedParameterStatus(params2);

    if(ii == 1 || mod(ii-1,500) == 0)
        m_mu = zeros(size(mu_f));
        v_mu = zeros(size(mu_f));
        m_L = zeros(size(L_f));
        v_L = zeros(size(L_f));
    end


    mu_step = mean(dlp,2);
    m_mu = beta_1 * m_mu + (1-beta_1)*mu_step;
    v_mu = beta_2 * v_mu + (1-beta_2)*mu_step.^2;
    m_hat_mu = m_mu./(1-beta_1^ii);
    v_hat_mu = v_mu./(1-beta_2^ii);
    mu_f = mu_f + alpha*m_hat_mu./(sqrt(v_hat_mu) + epsilon);
    
    L_step  = (dlp*eta')./M + diag(1./diag(L_f));
    m_L = beta_1 * m_L + (1-beta_1)*L_step;
    v_L = beta_2 * v_L + (1-beta_2)*L_step.^2;
    m_hat_L = m_L./(1-beta_1^ii);
    v_hat_L = v_L./(1-beta_2^ii);
    L_f = L_f + tril(alpha*m_hat_L./(sqrt(v_hat_L) + epsilon));
%     L = diag(diag(L));
    for kk = 1:size(L_f,1)
        L_f(kk,kk) = sign(L_f(kk,kk))*max(1e-4, abs(L_f(kk,kk)));
    end
    L_f(bb) = 0;


    nl_mu_prev = -nllFunction(mu_f);
    nl_mu(ii) = nl_mu_prev;

    elp(ii) = mean(lp);
    ell(ii) = mean(ll);
    ell_prev = ell(ii);

    if(mod(ii,20) == 0 || ii == 1)
        %%
        figure(2);
        clf
        subplot(1,4,1)
        plot([nl_mu ell]);
        subplot(1,4,2)
        plot(mu_f)
        subplot(1,4,3)
        plot(diag(L_f))
        subplot(1,4,4)
        plot(L_f(lls));
        drawnow;
        pause(1);
    end
end
