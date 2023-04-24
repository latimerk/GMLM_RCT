function [] = printRegularizedParameterStatus(params)
H = double(params.H);

K = size(params.B,1);
if(K > 0)
    if(numel(H) > 2)
        w_mu_idx  = 1;
        w_sig_idx = 2;
        b_mu_idx = 2 + (1:K);
        b_sig_idx = (2 + K + 1):numel(H);
        
        w_mu = H(w_mu_idx);
        b_mu = H(b_mu_idx);
    else
        w_sig_idx = 1;
        b_sig_idx = 2;
    end
    log_w_sig = H(w_sig_idx);
    log_b_sig = H(b_sig_idx);
    if(numel(log_b_sig) > 1)
        fprintf("\t\tw_sig = %.1e, min(b_sig) = %.1e, max(b_sig) = %.1e, ", exp(log_w_sig) + 0e-3, exp(min(log_b_sig)), exp(max(log_b_sig)));
    else
        fprintf("\t\tw_sig = %.1e, b_sig = %.1e, ", exp(log_w_sig) + 0e-3, exp(log_b_sig));
    end
else
    w_mu_idx  = 1;
    w_sig_idx = 2;
    log_w_sig = H(w_sig_idx);
    fprintf("\t\tw_sig = %.1e, ", exp(log_w_sig) + 0e-3);
end

for ii = 1:numel(params.Groups)
%     fprintf("c_%d = %.2f, tau_%d = %.2f, max(|V_%d|) = %.2f", ii, exp(1/2*params.Groups(ii).H(1)), ii, exp(params.Groups(ii).H(2)), ii, max(sqrt(sum(params.Groups(ii).V.^2,2))));

    if(numel(params.Groups(ii).T) == 2)
        r = sqrt(sum(params.Groups(ii).V.^2,1)).*sqrt(sum(params.Groups(ii).T{1}.^2,1)).*sqrt(sum(params.Groups(ii).T{2}.^2,1));
    else
        r = sqrt(sum(params.Groups(ii).V.^2,1)).*sqrt(sum(params.Groups(ii).T{1}.^2,1));
    end
    rank_c = size(params.Groups(ii).V,2);

    if(numel(params.Groups(ii).H) == 1)
        fprintf("tau_%d = %.1e, r_%d = [%.1e, %.1e, %.1e, R=%d]", ii, exp(1/2*params.Groups(ii).H(1)), ii, min(r), median(r), max(r), rank_c);
    else
        fprintf("c_%d = %.1e, tau_%d = %.1e, r_%d = [%.1e, %.1e, %.1e, R=%d]", ii, exp(1/2*params.Groups(ii).H(1)), ii, exp(params.Groups(ii).H(2)), ii, min(r), median(r), max(r), rank_c);
    end
    if(ii < numel(params.Groups))
        fprintf(", ");
    end
end
fprintf("\n");