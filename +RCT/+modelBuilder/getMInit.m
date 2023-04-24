function [M_setup] = getMInit(params, modelSetup, numTrials, M_weights)

M_setup = params;

if(~isfield(modelSetup, "useRidge") || ~modelSetup.useRidge)
    if(nargin < 4 || isempty(M_weights))
        M_weights.cs   = 1;
        M_weights.taus = 2;
        M_weights.phis = 1;
        M_weights.lambda_latent = 1/min(8,numTrials);
        M_weights.lambdas = 1;

        M_weights.lambdas_V = 1;
        M_weights.V = 1;%M_weights.lambdas_V/2;
    
        M_weights.W = 2;
        M_weights.B = 2;
        M_weights.H_W_mu = 2;
        M_weights.H_B_mu = 2;
        M_weights.H_W_sig = 2;
        M_weights.H_B_sig = 2;

        M_weights.Kernels = 1;
    end
    
    B = size(params.B,1);
    P = size(params.B,2);
    
    M_setup.W(:) = M_weights.W;

    M_setup.B(:,:) = M_weights.B;

    M_setup.H(1)   = M_weights.H_W_mu;
    M_setup.H(2)   = M_weights.H_W_sig;
    M_setup.H(2 + (1:B)) = M_weights.H_B_mu;
    M_setup.H((2 + B + 1):end) = M_weights.H_B_sig;
    
    
    jj = 1;
    R = size(params.Groups(jj).V,2);
    NV_lambdas = numel(params.Groups(jj).V);
    M_setup.Groups(jj).H(:) = M_weights.lambdas;
    M_setup.Groups(jj).H(1) = M_weights.cs;
    M_setup.Groups(jj).H(2) = M_weights.taus;
    M_setup.Groups(jj).H(2+(1:R)) = M_weights.phis;
    M_setup.Groups(jj).H(2+R+(1:NV_lambdas)) = M_weights.lambdas_V;

    M_setup.Groups(jj).V(:) = M_weights.V ;
    M_setup.Groups(jj).T{1}(:) = M_weights.Kernels;
    M_setup.Groups(jj).T{2}(1:size(modelSetup.stimConfig(1).regressors,2), :) = 1;
    M_setup.Groups(jj).T{2}((size(modelSetup.stimConfig(1).regressors,2)+1):end,:) = M_weights.lambda_latent;
    
    jj = 2;
    R = size(params.Groups(jj).V,2);
    NV_lambdas = numel(params.Groups(jj).V);
    M_setup.Groups(jj).H(:) = M_weights.lambdas;
    M_setup.Groups(jj).H(1) = M_weights.cs;
    M_setup.Groups(jj).H(2) = M_weights.taus;
    M_setup.Groups(jj).H(2+(1:R)) = M_weights.phis;
    M_setup.Groups(jj).V(:) = M_weights.V ;
    M_setup.Groups(jj).H(2+R+(1:NV_lambdas)) = M_weights.lambdas_V;
    M_setup.Groups(jj).T{1}(:) = M_weights.Kernels;

    M_setup.Groups(jj).T{2}(1:size(modelSetup.responseConfig(1).regressors,2), :) = 1;
    M_setup.Groups(jj).T{2}((size(modelSetup.responseConfig(1).regressors,2)+1):end,:) = M_weights.lambda_latent;
    
    for jj = 3:numel(params.Groups)
        R = size(params.Groups(jj).V,2);
        NV_lambdas = numel(params.Groups(jj).V);
        M_setup.Groups(jj).H(:) = M_weights.lambdas;
        M_setup.Groups(jj).H(1) = M_weights.cs;
        M_setup.Groups(jj).H(2) = M_weights.taus;
        M_setup.Groups(jj).H(2+(1:R)) = M_weights.phis;
        M_setup.Groups(jj).H(2+R+(1:NV_lambdas)) = M_weights.lambdas_V;
        M_setup.Groups(jj).V(:) = M_weights.V ;
        M_setup.Groups(jj).T{1}(:) = M_weights.Kernels;
        M_setup.Groups(jj).T{2}(:) = 1;

    end
else
    if(nargin < 4 || isempty(M_weights))
        M_weights.taus = 4;
    
        M_weights.W = 8;
        M_weights.B = 4;
        M_weights.H = 4;
    end
    
    
    M_setup.W(:) = M_weights.W;
    M_setup.B(:) = M_weights.B;
    M_setup.H(:)   = M_weights.H;
    
    jj = 1;
    M_setup.Groups(jj).H(:) = M_weights.taus;
    M_setup.Groups(jj).V(:) = 1;
    M_setup.Groups(jj).T{1}(:) = 1;
    M_setup.Groups(jj).T{2}(1:size(modelSetup.stimConfig(1).regressors,2), :) = 1;
    
    jj = 2;
    M_setup.Groups(jj).H(:) = M_weights.taus;
    M_setup.Groups(jj).V(:) = 1;
    M_setup.Groups(jj).T{1}(:) = 1;
    M_setup.Groups(jj).T{2}(1:size(modelSetup.responseConfig(1).regressors,2), :) = 1;
    
    for jj = 3:numel(params.Groups)
        M_setup.Groups(jj).H(:) = M_weights.taus;
        M_setup.Groups(jj).V(:) = 1;
        M_setup.Groups(jj).T{1}(:)    = 1;
        M_setup.Groups(jj).T{2}(:, :) = 1;
    end
end