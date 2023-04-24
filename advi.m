NI = 4000;
alphaSteps.step  = [1     1000  3500];
alphaSteps.alpha = [1e-5  1e-4   1e-4];
alphaSteps.reset = [1];
M = 50;

% ADVI/ BBVI: full rank

opts = gmlm.getComputeOptionsStruct(true);
res = gmlm.getEmptyResultsStruct(opts);
w_init = gmlm.vectorizeParams(params, opts);
P = numel(w_init);
nlpostFunction = @(ww) gmlm.vectorizedNLPost_func(ww, params, opts, res);
nllFunction = @(ww) gmlm.vectorizedNLL_func(ww, params, opts, res);

% mu_f = w_init;
% L_f  = 0.1*eye(P);
mu_f = mu;
L_f = diag(exp(L));

elp = nan(NI, 1);
ell = nan(NI, 1);
H = nan(NI, 1);
nl_mu = nan(NI, 1);
nl_mu_prev = -nllFunction(mu_f);
ell_prev = nl_mu_prev;
alpha = 0.001;
beta_1 = 0.9;
beta_2 = 0.999;
epsilon = 1e-8;


for ii = 1:NI
    fprintf("ii = %d / %d (%e, %e)\n", ii, NI,nl_mu_prev,ell_prev);
    eta  = randn(P, M);
    zeta = L_f*eta + mu_f;

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
    L_step  = tril(dlp*eta')./M + diag(1./diag(L_f));


%     tau = 1;
%     beta = 0.1; % eta in paper
%     alpha = 0.1;
%     epsilon = 1e-5;
%     if(ii == 1)
%         s_mu = mu_step.^2;
%         s_L  = L_step.^2;
%     else
%         s_mu = alpha*mu_step.^2 + (1-alpha).*s_mu;
%         s_L  = alpha*L_step.^2  + (1-alpha).*s_L;
%     end
%     p_mu = beta*ii^(-1/2+epsilon)./(tau + sqrt(s_mu));
%     p_L  = tril(beta*ii^(-1/2+epsilon)./(tau + sqrt(s_L)));
%     mu = mu + p_mu.*mu_step;
%     L  = L + p_L.*L_step;



    if(ismember(ii, alphaSteps.reset))
        m_mu = zeros(size(mu_f));
        v_mu = zeros(size(mu_f));
        m_L = tril(zeros(size(L_f)));
        v_L = tril(zeros(size(L_f)));
    end
    if(ismember(ii, alphaSteps.step))
        [~,jj] = find(alphaSteps.step == ii,1,"first");
        a_c = alphaSteps.alpha(jj);
        alpha_mu = a_c;
        alpha_L = a_c;
    end
    m_mu = beta_1 * m_mu + (1-beta_1)*mu_step;
    v_mu = beta_2 * v_mu + (1-beta_2)*mu_step.^2;
    m_hat_mu = m_mu./(1-beta_1^ii);
    v_hat_mu = v_mu./(1-beta_2^ii);
    
    m_L = (beta_1 * m_L + (1-beta_1)*L_step);
    v_L = (beta_2 * v_L + (1-beta_2)*L_step.^2);
    m_hat_L = (m_L./(1-beta_1^ii));
    v_hat_L = (v_L./(1-beta_2^ii));
    
    mu_f = mu_f + alpha_mu*m_hat_mu./(sqrt(v_hat_mu) + epsilon);
    L_f = L_f + alpha_L*(m_hat_L./(sqrt(v_hat_L) + epsilon));
    if(~isreal(L_f))
        error("WHAT");
    end

    elp(ii) = mean(lp);
    ell(ii) = mean(ll);
    H(ii) = sum(log(abs(diag(L_f))));
    nl_mu_prev = -nllFunction(mu_f);
    nl_mu(ii) = nl_mu_prev;
    ell_prev = ell(ii);

    if(mod(ii,50) == 0)
        kgmlm.utils.sfigure(1);
        clf
        subplot(1,3,1)
        if(ii > 1000) 
            plot(1000:NI, [nl_mu(1000:end) ell(1000:end)]);
        elseif(ii > 500) 
            plot(500:NI, [nl_mu(500:end) ell(500:end)]);
        elseif(ii > 200) 
            plot(200:NI, [nl_mu(200:end) ell(200:end)]);
        else
            plot([nl_mu ell]);
        end
        subplot(1,3,2)
        plot(mu_f)
        subplot(1,3,3)
        imagesc(L_f)
        drawnow;
        pause(0.1);
    end
end
%%

params_advi_0 = gmlm.devectorizeParams(mu_f, params, opts);
params_advi_r = gmlm.rescaleParamStruct(params_advi_0, scaled_WB, scaled_VT);

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
zeta = L_f*eta + mu_f;

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

save("adviRun_f.mat", "-v7.3", "params_advi_0", "params_advi_r", "mu_f", "L_f", "samples_advi_r");