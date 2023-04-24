function [prior_normalization, spkHist] = computeNormalizationSettings(trials, bases, kernelTimescales, t_pre_noise_bins, t_post_saccade_bins, couplingShuffle, baselineRates, baselineCorrectionType)
if(nargin < 6)
    couplingShuffle = [];
end
if(nargin < 8)
    baselineRates = [];
    baselineCorrectionType = "none";
end
if(baselineCorrectionType ~= "none" && isempty(baselineRates))
    error("Cannot correct baseline rates for spike history - no baseline rates given.");
end


NB_h = size(bases.spkHist.B,2);


regions = fieldnames(trials(1).Y);

regions_struct = cell(2, numel(regions));
for ff = 1:numel(regions)
    regions_struct{1,ff} = regions{ff};
end
spkHist = struct(regions_struct{:});
% spkTimes = struct(regions_struct{:});

if(NB_h > 0)
    B = bases.spkHist.B;
    paddedHspk     = [zeros(size(B)     + [1 0]); B];%zero padding is to make convolutions easier
    L = kernelTimescales.spkHist;
    tts = bases.spkHist.tts_0;
    stdParams = bases.spkHist.std;
    D = kernelFunction(tts, L, stdParams);
    S = B'*(D*B);
    
    prior_normalization.spkHist.U_sigma_0   = S;
    prior_normalization.spkHist.U_sigma     = eye(size(S));
    prior_normalization.spkHist.U_transform = chol(S);
    
    %% convolves spike history for each trial
    
    for tt = 1:numel(trials)
    
        t_0 = trials(tt).noise_on - t_pre_noise_bins;
        tts_c = (t_0):(trials(tt).saccade_end + t_post_saccade_bins);
        T_c = numel(tts_c);
    
        %% convolves spike history for each area
        for ff = 1:numel(regions)
            region_c = regions{ff};
    
            tts_c_coupling = tts_c;
            tt_c = tt;
            if(~isempty(couplingShuffle) && isfield(couplingShuffle, region_c) && ~isempty(couplingShuffle.(region_c)))
                tt_c = couplingShuffle.(region_c)(tt);
                if(tt_c ~= tt)
                    dt_0_coupling = trials(tt_c).(couplingShuffle.event.(region_c) ) - trials(tt).(couplingShuffle.event.(region_c) );
                    tts_c_coupling = tts_c + dt_0_coupling;
                    if(tts_c_coupling(1) < 1 || tts_c_coupling(end) > size(trials(tt_c).Y.(region_c),1))
                        error("Trial window exceeds processed data window size!");
                    end
                end
            end
    
            Y_t = trials(tt_c).Y.(region_c);
            
            if(baselineCorrectionType ~= "none")
                Y_b = baselineRates.(region_c).meanRates(tt_c).Y;
                Y_t = correctSpikeHistory(Y_t, Y_b, baselineCorrectionType);
            end

            dim_P_c = size(Y_t,2);
            spkHist(tt).(region_c) = zeros(T_c, NB_h, dim_P_c);
    
            for bb = 1:NB_h
                h_c = conv2(Y_t, paddedHspk(:, bb), "same");
                spkHist(tt).(region_c)(:, bb, :) = reshape(h_c(tts_c_coupling, :), T_c,1,[]);
            end
    
%             [spTimes,neuronNums] = find(Y_t);
%             spTimes = spTimes - tts_c_coupling(1) + 1;
%             spkTimes(tt).(region_c)  = cell(dim_P_c,1); 
%             for pp = 1:dim_P_c
%                 spkTimes(tt).(region_c){pp} = sort(spTimes(neuronNums == pp));
%             end
        end
    end

    %% gets mean & std info about convolved spike trains per neuron
    cS = prior_normalization.spkHist.U_transform;
    for ff = 1:numel(regions)
        region_c = regions{ff};
    
        xx = cell2mat(arrayfun(@(aa)aa.(region_c), spkHist', 'UniformOutput', false));
    
        std_per_neuron = nan(size(xx,3),1);
        for pp = 1:size(xx,3)
            std_per_neuron(pp) = sqrt(trace(cov(xx(:,:,pp) * cS)));% * sqrt(2/pi); %% xx = convolved spike history
        end
        prior_normalization.spkHist.(region_c).mu = mean(xx,1);
        prior_normalization.spkHist.(region_c).sig = std_per_neuron;
    end

    %% subtract means from colved spike history
    for ff = 1:numel(regions)
        region_c = regions{ff};
        for tt = 1:numel(trials)
            spkHist(tt).(region_c) = spkHist(tt).(region_c) - prior_normalization.spkHist.(region_c).mu;
        end
    end

end

%% per task basis
basisNames = fieldnames(bases);
for bb = 1:numel(basisNames)
    basis = basisNames{bb};
    if(strcmpi(basis, "spkHist") || ~isstruct(bases.(basis)))
        continue;
    end
    B = bases.(basis).B;
    L = kernelTimescales.(basis);
    tts = bases.(basis).tts_0;
    stdParams = bases.(basis).std;
    D = kernelFunction(tts, L, stdParams);
    S = B'*(D*B);
    
    prior_normalization.(basis).U_sigma_0   = S;
    prior_normalization.(basis).U_sigma     = eye(size(S));
    prior_normalization.(basis).U_transform = chol(S);
    
    %cS = prior_normalization.(basis).U_transform;
    std_per_basis  = 1;%sqrt(trace(cov(B * cS))) * sqrt(2/pi); % not sure why I was computing this before
    
    prior_normalization.(basis).mu = 0;
    prior_normalization.(basis).sig = std_per_basis;
end

end

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


function [Y_c] = correctSpikeHistory(Y_t, Y_b, correctionType)
    switch(correctionType)
        case "baseline"
            Y_c = Y_b;
        case "subtract"
            Y_c =  Y_t - Y_b;
        case "rate"
            Y_c =  Y_b;
        case "zscore"
            Y_c =  (Y_t - Y_b)./sqrt(Y_b.*(1-Y_b));
        case "bernoulli"
            Y_c =  double(rand(size(Y_b)) < Y_b);
        case "none"
            Y_c = Y_t;
        otherwise
            error("Unknown baseline correction time");
    end
    Y_c(isnan(Y_c)) = 0;
end