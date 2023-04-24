addpath ../GMLM_RCT/
addpath ../GMLM
addpath ../GMLM/example/
bases = RCT.modelBuilder.setupBasis();
kernelTimescales =  RCT.modelBuilder.getDefaultKernelTimescales();

basis = "stimulus";

B = bases.(basis).B;
L = kernelTimescales.(basis);
tts = bases.(basis).tts_0;
stdParams = bases.(basis).std;
D = kernelFunction(tts, L, stdParams);
S = B'*(D*B);

prior_normalization.(basis).U_sigma_0   = S;
prior_normalization.(basis).U_sigma     = eye(size(S));
prior_normalization.(basis).U_transform = chol(S);

cS = prior_normalization.(basis).U_transform;
std_per_basis  = sqrt(trace(cov(B * cS))) * sqrt(2/pi);

prior_normalization.(basis).mu = 0;
prior_normalization.(basis).sig = std_per_basis;

function [C] =  kernelFunction(tts, L, stdSettings)
D = tts(:) - tts(:)';

% matern kernel (3/2)
C_0 = (1 + sqrt(3)*abs(D)./L ).*exp(-sqrt(3)*abs(D)./L);

% matern kernel (5/2)
% C_0 = (1 + sqrt(5)*abs(D)./L + 5*D.^2./(3*L^2)).*exp(-sqrt(5)*abs(D)./L);

% OU
% C_0 = exp(-abs(D)./L);
% White noise
% C_0 = eye(size(D,1));

if(nargin > 2 && ~isempty(stdSettings))
    S = RCT.modelBuilder.getSTDForKernels(tts, stdSettings);
    C = C_0.*(S(:)*S(:)');

    if(~all(stdSettings.mask))
        C(~stdSettings.mask, :) = 0;
        C(:, ~stdSettings.mask) = 0;
        C(~stdSettings.mask, ~stdSettings.mask) = eye(sum(~stdSettings.mask)) * stdSettings.mask_std.^2;
    end
else
    C = C_0;
end

end