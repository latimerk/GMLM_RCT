function [GMLMstructure, trials, spkHistPrior_setup, prior_setup, trialsUsed, couplingShuffle, prior_normalization, timeWindows] = constructGMLMdata14(data, modelSetup, varargin)

p = inputParser;
p.CaseSensitive = false;

addRequired(p, "modelSetup"   ,    @(A) isstruct(A) );

default_noise_time = 0;
includeFixation = isfield(modelSetup, "fixationConfig") && ~isempty(modelSetup.fixationConfig);
if(includeFixation)
    default_noise_time = 400;
end

addParameter(p, "t_pre_noise_bins"   ,  default_noise_time  / data.bin_size_ms,    @isnumeric);
addParameter(p, "t_post_saccade_bins",  60 / data.bin_size_ms,    @isnumeric);
addParameter(p, "baselineRates",  [],    @(aa)isempty(aa) || isstruct(aa));

validate_data = @(A) isstruct(A) && isfield(A, "NeuronInfo") && isfield(A,"TaskInfo") && isfield(A, "trials") && isstruct(A.trials) && isfield(A.trials,"Y");
if(~validate_data(data))
    error("Invalid data struct");
end

parse(p, modelSetup, varargin{:});
% then set/get all the inputs out of this structure
t_pre_noise_bins  = p.Results.t_pre_noise_bins;
t_post_saccade_bins  = p.Results.t_post_saccade_bins;
couplingShuffle = [];


if(~isfield(modelSetup, "runNumber"))
    modelSetup.runNumber = 0;
end

baselineRates = p.Results.baselineRates;
if(~isfield(modelSetup, "baselineCorrectionType") || isempty(baselineRates))
    modelSetup.baselineCorrectionType = "none";
end

%% check if requested recording location exists
modelSetup.location = string(modelSetup.location).upper;
if(~isfield(data.trials(1).Y, modelSetup.location))
    error("Recording location '%s' not found!", modelSetup.location);
end

locations = ["FEF", "LIP", "SC"];
[~,loc_main_idx] = ismember(modelSetup.location, locations);
if(isempty(loc_main_idx))
    error("Unknown coupling location");
end

if(modelSetup.targets_verical)
    target_str = "vertical target";
else
    target_str = "horizontal target";
end
correctly_oriented_trials = [data.trials(:).targets_vertical] == modelSetup.targets_verical;


trialsUsed = find(correctly_oriented_trials);
trialsUsed = trialsUsed(:);

qtr1 = trialsUsed(1:4:numel(trialsUsed));
qtr2 = trialsUsed(2:4:numel(trialsUsed));
qtr3 = trialsUsed(3:4:numel(trialsUsed));
qtr4 = trialsUsed(4:4:numel(trialsUsed));

half1 = trialsUsed(1:2:numel(trialsUsed));
half2 = trialsUsed(2:2:numel(trialsUsed));

if(isfield(modelSetup, "fitQuarter") && all(modelSetup.fitQuarter > 0))
    trialsUsed = [];
    if(ismember(1,modelSetup.fitQuarter))
        trialsUsed = [trialsUsed;qtr1];
    end
    if(ismember(2,modelSetup.fitQuarter))
        trialsUsed = [trialsUsed;qtr2];
    end
    if(ismember(3,modelSetup.fitQuarter))
        trialsUsed = [trialsUsed;qtr3];
    end
    if(ismember(4,modelSetup.fitQuarter))
        trialsUsed = [trialsUsed;qtr4];
    end
    if(isempty(trialsUsed))
        error("invalid fitQuarter setting.");
    end
elseif(isfield(modelSetup, "fitHalf") && modelSetup.fitHalf > 0)
    if(modelSetup.fitHalf == 1)
        trialsUsed = half1;
    elseif(modelSetup.fitHalf == 2)
        trialsUsed = half2;
    else
        error("invalid fitHalf setting.");
    end
end
trialsUsed = sort(trialsUsed);
data.trials = data.trials(trialsUsed);
dim_M = numel(data.trials);

if(modelSetup.baselineCorrectionType ~= "none" && ~isempty(baselineRates))
    for ii = 1:numel(data.NeuronInfo)
        loc_c = data.NeuronInfo(ii).location;
        baselineRates.(loc_c).meanRates = baselineRates.(loc_c).meanRates(trialsUsed);
    end
end

%% do any trial shuffling
if(isfield(modelSetup, "interAreaCoupling_shuffle") && modelSetup.interAreaCoupling_shuffle)
    for ii = 1:numel(data.NeuronInfo)
        loc_c = data.NeuronInfo(ii).location;
    
        if(~strcmpi(loc_c,modelSetup.location) && isfield(modelSetup, "interAreaCoupling_shuffle") && modelSetup.interAreaCoupling_shuffle)
            rng_0 = rng();
            [~,loc_c_idx] = ismember(loc_c, locations);
            if(isempty(loc_c_idx))
                error("Unknown coupling location");
            end
            rng(modelSetup.runNumber + 1000 * loc_main_idx + 100*loc_c_idx);
            couplingShuffle.(loc_c) = randperm(dim_M);
            rng(rng_0);
    
            if(modelSetup.interAreaCoupling_shuffle == 2)
                couplingShuffle.event.(loc_c) = "saccade_end";
            else
                couplingShuffle.event.(loc_c) = "noise_on";
            end
        else
            couplingShuffle.(loc_c) = [];
            couplingShuffle.event.(loc_c) = "";
        end
    end
else
    couplingShuffle = [];
end

%%
[prior_normalization, spkHist] = RCT.modelBuilder.computeNormalizationSettings(data.trials, modelSetup.bases, modelSetup.kernelTimescales, t_pre_noise_bins, t_post_saccade_bins, couplingShuffle, baselineRates, modelSetup.baselineCorrectionType);
if(isfield(modelSetup, "normalizationSettings") && ~isempty(modelSetup.normalizationSettings))
    prior_normalization = modelSetup.normalizationSettings;
end

%% Setup GMLM groups
GMLMstructure.dim_B = size(modelSetup.bases.spkHist.B, 2); % specify size of the linear term
GMLMstructure.dim_P = size(data.trials(1).Y.(modelSetup.location), 2);
GMLMstructure.Groups = struct("X_shared", [], "dim_R_max", [], "dim_A", [], "name", [], "dim_names", [], "dim_T", [], "dim_F", [], "factor_idx", [], "is_coupling", []); % tensor coefficient groups

%  stimulus setup (DENSE / local regressors)

eventsStimulus = ["noise_on";
                    "stim_on";
                    "targets_on";
                    "fix_off"];
eventsStimulus = eventsStimulus(ismember(eventsStimulus, [modelSetup.stimConfig(:).event_name]));

eventsResponse = "saccade_end";

[~,event_idx] =  ismember(eventsStimulus, [modelSetup.stimConfig(:).event_name]);
%[~,event_idx_r] = ismember(eventsResponse, [modelSetup.responseConfig(:).event_name]);
if(includeFixation)
    eventsFixation = "noise_on";
    %[~,event_idx_f] = ismember(eventsFixation, [modelSetup.fixationConfig(:).event_name]);
end
BB = cell2mat({modelSetup.stimConfig(event_idx).regressors}');
cond_init = cumsum([0;arrayfun(@(a) size(a.regressors,1), modelSetup.stimConfig(event_idx) )]);

localCouplingGroup = [];

if(~modelSetup.includeTrialwiseLatent)
    jj = 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Stimulus";
    GMLMstructure.Groups(jj).dim_names = ["timing", "xstim"]; %dimensions of the tensor
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.stimulus; % max allocated space for rank
    GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.stimulus.B, 2) size(modelSetup.stimConfig(1).regressors,2)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    
    GMLMstructure.Groups(jj).dim_A = numel(eventsStimulus);    
    GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.stimulus.B, BB};
    GMLMstructure.Groups(jj).dim_F = (GMLMstructure.Groups(jj).dim_T);
    GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
    GMLMstructure.Groups(jj).grp_coef_names =  [modelSetup.stimConfig(:).coeff_names];
    GMLMstructure.Groups(jj).grp_coef_idxs  =  [modelSetup.stimConfig(:).coeff_groups];
    GMLMstructure.Groups(jj).grp_log_constant_scales =  zeros(1, numel(GMLMstructure.Groups(jj).grp_coef_idxs));
    GMLMstructure.Groups(jj).loadings_log_constant_scales = 0;
    GMLMstructure.Groups(jj).extra_log_constant_scale = 0;

    % Response setup
    jj = jj + 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Response";
    GMLMstructure.Groups(jj).dim_A = 1;     
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.response; % max allocated space for rank
    if(size(modelSetup.responseConfig.regressors,2) > 1) % if a tensor is necessary
        GMLMstructure.Groups(jj).dim_names = ["timing", "xchoice"]; %dimensions of the tensor
        GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.response.B, 2) size(modelSetup.responseConfig(1).regressors,2)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    
        GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.response.B, modelSetup.responseConfig(1).regressors};
        GMLMstructure.Groups(jj).dim_F = GMLMstructure.Groups(jj).dim_T;
        GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
        GMLMstructure.Groups(jj).grp_coef_names =  [modelSetup.responseConfig(:).coeff_names];
        GMLMstructure.Groups(jj).grp_coef_idxs  =  [modelSetup.responseConfig(:).coeff_groups];
        GMLMstructure.Groups(jj).grp_log_constant_scales =  zeros(1, numel(GMLMstructure.Groups(jj).grp_coef_idxs));
    
    else %low rank matrix instead
        GMLMstructure.Groups(jj).dim_names = "timing"; %dimensions of the tensor
        GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.response.B };
        GMLMstructure.Groups(jj).dim_T = size(modelSetup.bases.response.B, 2) ; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
        GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
        GMLMstructure.Groups(jj).factor_idx = 1; %factor setup
        GMLMstructure.Groups(jj).grp_coef_names =  [];
        GMLMstructure.Groups(jj).grp_coef_idxs  =  [];
        GMLMstructure.Groups(jj).grp_log_constant_scales = [];
    end

    GMLMstructure.Groups(jj).loadings_log_constant_scales = 0;
    GMLMstructure.Groups(jj).extra_log_constant_scale = 0;

    % fixation setup
    if(includeFixation)
        fprintf("Including fixation event.\n");
        jj = jj + 1;
        GMLMstructure.Groups(jj).is_coupling = false;
        GMLMstructure.Groups(jj).name = "Fixation";
        GMLMstructure.Groups(jj).dim_A = 1;     
        GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.fixation; % max allocated space for rank
        if(size(modelSetup.fixationConfig.regressors,2) > 1) % if a tensor is necessary
            GMLMstructure.Groups(jj).dim_names = ["timing", "taskDetails"]; %dimensions of the tensor
            GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.fixation.B, 2) size(modelSetup.fixationConfig(1).regressors,2)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
        
            GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.fixation.B, modelSetup.fixationConfig(1).regressors};
            GMLMstructure.Groups(jj).dim_F = GMLMstructure.Groups(jj).dim_T;
            GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
            GMLMstructure.Groups(jj).grp_coef_names =  [modelSetup.responseConfig(:).coeff_names];
            GMLMstructure.Groups(jj).grp_coef_idxs  =  [modelSetup.responseConfig(:).coeff_groups];
            GMLMstructure.Groups(jj).grp_log_constant_scales =  zeros(1, numel(GMLMstructure.Groups(jj).grp_coef_idxs));
        
        else %low rank matrix instead
            GMLMstructure.Groups(jj).dim_names = "timing"; %dimensions of the tensor
            GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.fixation.B };
            GMLMstructure.Groups(jj).dim_T = size(modelSetup.bases.fixation.B, 2) ; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
            GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
            GMLMstructure.Groups(jj).factor_idx = 1; %factor setup
            GMLMstructure.Groups(jj).grp_coef_names =  [];
            GMLMstructure.Groups(jj).grp_coef_idxs  =  [];
            GMLMstructure.Groups(jj).grp_log_constant_scales = [];
        end

        GMLMstructure.Groups(jj).loadings_log_constant_scales = 0;
        GMLMstructure.Groups(jj).extra_log_constant_scale = 0;
    end
else
    latentScaling = 1;

    if(modelSetup.includeTrialwiseLatent_allStimEvents)
        stim_latent_groups       = arrayfun(@(aa) size(aa.regressors,1), modelSetup.stimConfig(event_idx));
        stim_latent_groups_names = [modelSetup.stimConfig(event_idx).event_name];
    else
        stim_latent_groups = sum(arrayfun(@(aa) size(aa.regressors,1), modelSetup.stimConfig(event_idx)));
        stim_latent_groups_names = "stim";
    end

    jj = 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Stimulus";
    GMLMstructure.Groups(jj).dim_names = ["timing", "xstim_trial"]; %dimensions of the tensor
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.stimulus; % max allocated space for rank
    GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.stimulus.B, 2) size(BB,2)+dim_M*numel(stim_latent_groups)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct

    MM = zeros(sum(stim_latent_groups) * dim_M, dim_M * numel(stim_latent_groups));
    for ii = 1:numel(stim_latent_groups)
        bb = zeros(sum(stim_latent_groups),1);
        bb((cond_init(ii)+1):(cond_init(ii+1))) = 1;

        MM(:, (1:dim_M) + (ii-1)*dim_M) = kron(latentScaling*eye(dim_M), bb);
    end
    
    GMLMstructure.Groups(jj).dim_A = numel(eventsStimulus);    
    GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.stimulus.B, [repmat(BB, dim_M, 1) MM]};
    GMLMstructure.Groups(jj).dim_F = (GMLMstructure.Groups(jj).dim_T);
    GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup

    cn = [modelSetup.stimConfig(:).coeff_names];
    ci = [modelSetup.stimConfig(:).coeff_groups];
    GMLMstructure.Groups(jj).grp_log_constant_scales =  [zeros(1, numel(cn)) repmat(modelSetup.includeTrialwiseLatent_log_scaleFunc(dim_M), 1, numel(stim_latent_groups))];
    for ii = 1:numel(stim_latent_groups)
        cn = [cn sprintf("trial_%s", stim_latent_groups_names(ii))]; %#ok<AGROW> 
        ci = [ci (size(BB,2)+ (ii-1)*dim_M + (1:dim_M) )]; %#ok<AGROW> 
    end
    GMLMstructure.Groups(jj).grp_coef_names =  cn;
    GMLMstructure.Groups(jj).grp_coef_idxs  =  ci;

    GMLMstructure.Groups(jj).loadings_log_constant_scales = 0;
    GMLMstructure.Groups(jj).extra_log_constant_scale = 0;
    
    % Response setup
    jj = jj + 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Response";
    GMLMstructure.Groups(jj).dim_A = 1;     
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.response; % max allocated space for rank

    GMLMstructure.Groups(jj).dim_names = ["timing", "xchoice_trial"]; %dimensions of the tensor
    BB_r = modelSetup.responseConfig(1).regressors;
    GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.response.B, 2) size(BB_r,2)+dim_M]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct

    GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.response.B, [repmat(BB_r, dim_M, 1) kron(latentScaling * eye(dim_M), ones(size(BB_r,1),1))]};
    GMLMstructure.Groups(jj).dim_F = GMLMstructure.Groups(jj).dim_T;
    GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup

    cn = [modelSetup.responseConfig(:).coeff_names];
    ci = [modelSetup.responseConfig(:).coeff_groups];
    GMLMstructure.Groups(jj).grp_log_constant_scales =  [zeros(1, numel(cn)) modelSetup.includeTrialwiseLatent_log_scaleFunc(dim_M)];
    cn = [cn "trials_1"];
    ci = [ci (size(BB_r,2) + (1:dim_M))];
    GMLMstructure.Groups(jj).grp_coef_names =  cn;
    GMLMstructure.Groups(jj).grp_coef_idxs  =  ci;

    GMLMstructure.Groups(jj).loadings_log_constant_scales = 0;
    GMLMstructure.Groups(jj).extra_log_constant_scale = 0;

    % fixation setup
    if(includeFixation)
        fprintf("Including fixation event.\n");
        jj = jj + 1;
        GMLMstructure.Groups(jj).is_coupling = false;
        GMLMstructure.Groups(jj).name = "Fixation";
        GMLMstructure.Groups(jj).dim_A = 1;     
        GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.fixation; % max allocated space for rank
    
        GMLMstructure.Groups(jj).dim_names = ["timing", "taskDetails"]; %dimensions of the tensor
        BB_f = modelSetup.fixationConfig(1).regressors;
        GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.fixation.B, 2) size(BB_f,2)+dim_M]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    
        GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.fixation.B, [repmat(BB_f, dim_M, 1) kron(latentScaling * eye(dim_M), ones(size(BB_f,1),1))]};
        GMLMstructure.Groups(jj).dim_F = GMLMstructure.Groups(jj).dim_T;
        GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup

        cn = [modelSetup.fixationConfig(:).coeff_names];
        ci = [modelSetup.fixationConfig(:).coeff_groups];
        GMLMstructure.Groups(jj).grp_log_constant_scales =  [zeros(1, numel(cn)) modelSetup.includeTrialwiseLatent_log_scaleFunc(dim_M)];
        cn = [cn "trials_1"];
        ci = [ci (size(BB_f,2) + (1:dim_M))];
        GMLMstructure.Groups(jj).grp_coef_names =  cn;
        GMLMstructure.Groups(jj).grp_coef_idxs  =  ci;

        GMLMstructure.Groups(jj).loadings_log_constant_scales = 0;
        GMLMstructure.Groups(jj).extra_log_constant_scale = 0;
    end
end

% Coupling
coupling_locations = [];
if(modelSetup.includeCouplingLocal)
    % local coupling
    jj = jj + 1;
    localCouplingGroup = jj;

    coupling_locations = cat(2, coupling_locations, modelSetup.location);
    GMLMstructure.Groups(jj).is_coupling = true;
    GMLMstructure.Groups(jj).coupling_location = modelSetup.location;
    GMLMstructure.Groups(jj).name = "Coupling_Local";
    GMLMstructure.Groups(jj).dim_names = ["kernel", "neuron_weights"]; %dimensions of the tensor
    GMLMstructure.Groups(jj).dim_A = 1;     
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.coupling_local; % max allocated space for rank
    GMLMstructure.Groups(jj).X_shared = {[]};
    GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.spkHist.B, 2) GMLMstructure.dim_P]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
    GMLMstructure.Groups(jj).factor_idx = [1 1]; %factor setup
    GMLMstructure.Groups(jj).grp_coef_names =  [];
    GMLMstructure.Groups(jj).grp_coef_idxs  =  [];

    GMLMstructure.Groups(jj).grp_log_constant_scales      =  -log(prior_normalization.spkHist.(modelSetup.location).sig);
    GMLMstructure.Groups(jj).loadings_log_constant_scales =  0;
    GMLMstructure.Groups(jj).extra_log_constant_scale     =  -0.5*log(GMLMstructure.Groups(jj).dim_T(2));% SCALE BY NUMBER OF NEURONS TOO! The square root is because this gets multiplied in twice: T{2} and V
end
if(modelSetup.includeCouplingInter)
    % from other areas
    for ii = 1:numel(data.NeuronInfo)
        loc_c = data.NeuronInfo(ii).location;
        if(~strcmpi(loc_c, modelSetup.location))



            jj = jj + 1;
            coupling_locations = cat(2, coupling_locations, loc_c);
            GMLMstructure.Groups(jj).is_coupling = true;
            GMLMstructure.Groups(jj).coupling_location = loc_c;
            dim_P_c = size(data.trials(1).Y.(loc_c),2);
            GMLMstructure.Groups(jj).name = sprintf("Coupling_%s", loc_c);
            GMLMstructure.Groups(jj).dim_names = ["kernel", "neuron_weights"]; %dimensions of the tensor
            GMLMstructure.Groups(jj).dim_A = 1;     
            GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.coupling_inter; % max allocated space for rank
            GMLMstructure.Groups(jj).X_shared = {[]};
            GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.spkHist.B, 2) dim_P_c]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
            GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
            GMLMstructure.Groups(jj).factor_idx = [1 1]; %factor setup
            GMLMstructure.Groups(jj).grp_coef_names =  [];
            GMLMstructure.Groups(jj).grp_coef_idxs  =  [];
            GMLMstructure.Groups(jj).grp_log_constant_scales      =  -log(prior_normalization.spkHist.(loc_c).sig); 
            GMLMstructure.Groups(jj).loadings_log_constant_scales =  0;
            GMLMstructure.Groups(jj).extra_log_constant_scale     =  -0.5*log(GMLMstructure.Groups(jj).dim_T(2));% SCALE BY NUMBER OF NEURONS TOO! The square root is because this gets multiplied in twice: T{2} and V
                
        end
    end
end
N_tensor_groups = numel(GMLMstructure.Groups);

%% sets up trials
trials = struct("Y", cell(numel(data.trials),1), "X_lin", [], "neuron", [], "Groups", []);

N_bins_response = 0;
N_bins_fixation = 0;
N_bins_stim     = 0;
N_bins_trial    = nan(numel(data.trials),1);
Spk_count_fixation = zeros(1, GMLMstructure.dim_P);
Spk_count_response = zeros(1, GMLMstructure.dim_P);
Spk_count_stim = zeros(1, GMLMstructure.dim_P);
Spk_count_trial = zeros(numel(data.trials), GMLMstructure.dim_P);

fprintf("Loading %d trials from %d cells in %s: %s\n", numel(data.trials), GMLMstructure.dim_P, modelSetup.location, target_str);
if(modelSetup.includeCouplingLocal)
    fprintf("\t Including local coupling.\n");
else
    fprintf("\t Local coupling OFF.\n");
end
if(modelSetup.includeCouplingInter && numel(coupling_locations) > 1)
    fprintf("\t Inter-area coupling: ");
    for ll = 2:numel(coupling_locations)
        fprintf("%s ", coupling_locations(ll));
    end
        fprintf("\n");
else
    fprintf("\t Inter-area coupling OFF.\n");
end

%%
timeWindows = cell(numel(data.trials),1);
for tt = 1:numel(data.trials)
    if(tt == 1 || mod(tt,50) == 0)
        fprintf("\t\t trial %d / %d...\n", tt, numel(data.trials));
    end
    t_0 = data.trials(tt).noise_on - t_pre_noise_bins;
    tts_c = (t_0):(data.trials(tt).saccade_end + t_post_saccade_bins);
    if(tts_c(1) < 1 || tts_c(end) > size(data.trials(tt).Y.(modelSetup.location),1))
        error("Trial window exceeds processed data window size!");
    end

    timeWindows{tt} = tts_c;

    % save spike count
    trials(tt).Y = data.trials(tt).Y.(modelSetup.location)(tts_c, :);
    TT = size(trials(tt).Y,1);
    N_bins_trial(tt) = TT;
    Spk_count_trial(tt,:) = sum(trials(tt).Y);
        
    % gets spike history
    if(GMLMstructure.dim_B > 0)
        trials(tt).X_lin = spkHist(tt).(modelSetup.location);
    else
        trials(tt).X_lin = [];
    end

    trials(tt).Groups = struct("X_local", cell(N_tensor_groups,1), "iX_shared", []);

    N_bins_stim = N_bins_stim + TT;
    Spk_count_stim = Spk_count_stim + sum(trials(tt).Y, 1);

    %% setup trial-wise fluctuations
    jj_0 = 2;
    
    %% setup coupling
    for cc = 1:numel(coupling_locations)
        jj = jj_0 + cc;
        loc_c = coupling_locations(cc);

        NF_c = numel(unique(GMLMstructure.Groups(jj).factor_idx));
        trials(tt).Groups(jj).X_local    = cell(1,NF_c); 
        trials(tt).Groups(jj).iX_shared  = cell(1,NF_c);

        trials(tt).Groups(jj).X_local{1} = spkHist(tt).(loc_c); 
    end

    %% setup stimulus
    correct_target_location = data.trials(tt).target_location; % location of CORRECT target
    if(data.trials(tt).targets_vertical)
        correct_target_location = correct_target_location(2);
    else
        correct_target_location = correct_target_location(1);
    end
    if(correct_target_location > 0)
        correct_target_location = 2; % right or up
    else
        correct_target_location = 1; % left or down
    end

    cat_c = data.TaskInfo.categories(data.trials(tt).direction_idx); % cat 1 = green, cat 2 = red
    if(cat_c == 1 && correct_target_location == 1)
        target_config = 2; % cat 1 (green) is left
    elseif(cat_c == 2 && correct_target_location == 1)
        target_config = 1; % cat 2 (red) is left
    elseif(cat_c == 1 && correct_target_location == 2)
        target_config = 1;
    elseif(cat_c == 2 && correct_target_location == 2)
        target_config = 2;
    end
    %rows: 1  = left/down  red,   up/right  green
    %      2  = left/down  green, up/right  red

    target_choice_location = data.trials(tt).saccade_location;
    if(data.trials(tt).targets_vertical)
        target_choice_location = target_choice_location(2);
    else
        target_choice_location = target_choice_location(1);
    end
    if(target_choice_location > 0)
        target_choice_location = 2; % right or up
    else
        target_choice_location = 1; % left or down
    end

    if(target_choice_location == 1 && target_config == 1)
        target_choice_color = 2;
    elseif(target_choice_location == 1 && target_config == 2)
        target_choice_color = 1;
    elseif(target_choice_location == 2 && target_config == 1)
        target_choice_color = 1;
    elseif(target_choice_location == 2 && target_config == 2)
        target_choice_color = 2;
    end

    %rows: 1  = left/down  red
    %      2  = left/down  green
    %      3  = right/up red
    %      4  = right/up green

    conditions = [1; %noise on
                  data.trials(tt).direction_idx; %stim on
                  target_config; %targets_on
                  target_config];%fix_off 
    
    jj = 1;  
    trials(tt).Groups(jj).iX_shared = {[]}; 
    iX  = zeros(TT, numel(eventsStimulus)); 
    iX2 = zeros(TT, numel(eventsStimulus)); 

    for aa = 1:numel(eventsStimulus)
        timing = data.trials(tt).(eventsStimulus(aa)) - t_0 + 1;
        if(isnan(timing))
            continue;
        end

        tts_event = timing + modelSetup.bases.stimulus.tts;
        vv = tts_event > 0 & tts_event <= TT;
        if(all(~vv))
            continue;
        end

        iX( tts_event(vv), aa) = modelSetup.bases.stimulus.tts_idx(vv); 
        iX2(tts_event(vv), aa) = conditions(aa) + cond_init(event_idx(aa));
    end
    if(~modelSetup.includeTrialwiseLatent)  
        trials(tt).Groups(jj).iX_shared = {iX, iX2};
        trials(tt).Groups(jj).X_local = {[], []};
    else  
        pp = ones(size(iX2)) * size(BB,1)*(tt-1);
        pp(iX2 <= 0) = 0;
        trials(tt).Groups(jj).iX_shared = {iX, iX2 + pp};
        trials(tt).Groups(jj).X_local = {[], []};
    end

    %% setup response
    conditions = target_choice_color + (target_choice_location-1)*2;

    jj = 2;

    iX  = zeros(TT, numel(eventsResponse)); 
    if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
        iX2 = zeros(TT, numel(eventsResponse)); 
    end

    for aa = 1:numel(eventsResponse)
        timing = data.trials(tt).(eventsResponse(aa)) - t_0 + 1;
        if(isnan(timing))
            continue;
        end

        tts_event = timing + modelSetup.bases.response.tts;
        vv = tts_event > 0 & tts_event <= TT;

        if(all(~vv))
            continue;
        end

        iX(tts_event(vv), aa) = modelSetup.bases.response.tts_idx(vv); 
        if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
            iX2(tts_event(vv), aa) = conditions(aa);
        end

        if(aa == 1)
            % get num bins and mean spike count contributing to response window
            N_bins_response = N_bins_response + numel(tts_event(vv));
            Spk_count_response = Spk_count_response + sum(trials(tt).Y(tts_event(vv),:),1);
        end
    end
    if(~modelSetup.includeTrialwiseLatent)  
        trials(tt).Groups(jj).iX_shared = {iX, iX2};
        trials(tt).Groups(jj).X_local = {[], []};
    else
        pp = ones(size(iX2)) * size(BB_r,1)*(tt-1);
        pp(iX2 <= 0) = 0;
        trials(tt).Groups(jj).iX_shared = {iX, iX2 + pp};
        trials(tt).Groups(jj).X_local = {[], []};
    end

    %% setup response
    if(includeFixation)
        conditions = 1;
        jj = 3;
    
        iX  = zeros(TT, numel(eventsFixation)); 
        if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
            iX2 = zeros(TT, numel(eventsFixation)); 
        end

        for aa = 1:numel(eventsFixation)
            timing = data.trials(tt).(eventsFixation(aa)) - t_0 + 1;
            if(isnan(timing))
                continue;
            end
    
            tts_event = timing + modelSetup.bases.fixation.tts;
            vv = tts_event > 0 & tts_event <= TT;
    
            if(all(~vv))
                continue;
            end
    
            iX(tts_event(vv), aa) = modelSetup.bases.fixation.tts_idx(vv); 
            if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
                iX2(tts_event(vv), aa) = conditions(aa);
            end
    
            if(aa == 1)
                % get num bins and mean spike count contributing to response window
                N_bins_fixation = N_bins_fixation + numel(tts_event(vv));
                Spk_count_fixation = Spk_count_fixation + sum(trials(tt).Y(tts_event(vv),:),1);
            end
        end
        if(~modelSetup.includeTrialwiseLatent)  
            if(numel(GMLMstructure.Groups(jj).dim_F) > 1)
                trials(tt).Groups(jj).iX_shared = {iX, iX2};
                trials(tt).Groups(jj).X_local = {[], []};
            else
                trials(tt).Groups(jj).iX_shared = {iX};
                trials(tt).Groups(jj).X_local = {[]};
            end
        else
            error("Not implemented");
        end
    end
end


%%
fprintf("Done.\n");

%% setup Horseshoe priors to encourage sparsity
spkHistPrior_setup.hyperprior.w_sig_nu = 1;
spkHistPrior_setup.hyperprior.b_sig_nu = 1;
spkHistPrior_setup.hyperprior.log_w_sig_scale = log(2);
spkHistPrior_setup.hyperprior.log_b_sig_scale = log(2);
spkHistPrior_setup.hyperprior.log_w_mu_sig = log(4);
spkHistPrior_setup.hyperprior.log_b_mu_sig = log(4);
spkHistPrior_setup.hyperprior.min_w_sig = 1e-2;
spkHistPrior_setup.hyperprior.min_b_sig = 1e-2;


if(GMLMstructure.dim_B > 0)
    spkHistPrior_setup.NH              = 2 + 1 + GMLMstructure.dim_B; %number of hyperparameters

    spkHistPrior_setup.U_sigma_0        = prior_normalization.spkHist.U_sigma_0;
    spkHistPrior_setup.U_sigma          = prior_normalization.spkHist.U_sigma;
    spkHistPrior_setup.U_transform      = prior_normalization.spkHist.U_transform;
    spkHistPrior_setup.U_sigma_inv      = inv(prior_normalization.spkHist.U_sigma_0);
    spkHistPrior_setup.U_sigma_chol_inv = chol(spkHistPrior_setup.U_sigma_inv)';
else
    spkHistPrior_setup.NH              = 2; %number of hyperparameters
end

if(~isfield(modelSetup, "flatPrior") || modelSetup.flatPrior)
    GMLMstructure.prior.log_prior_func = @(params, results, priorOnly) DMC.priors.Horseshoe3.spkHistPrior3(params, results, spkHistPrior_setup, priorOnly);
    GMLMstructure.scaleParams          = @(params) DMC.priors.Horseshoe3.scaleWB3(params, spkHistPrior_setup, localCouplingGroup);
    GMLMstructure.scaleDerivatives     = @(results,params,posterior,A) DMC.priors.Horseshoe3.scaleDWB3(results, params, spkHistPrior_setup, posterior, localCouplingGroup, A);
    GMLMstructure.prior.dim_H          = spkHistPrior_setup.NH;
else
    GMLMstructure.prior.log_prior_func = @(params, results, priorOnly) DMC.priors.HorseshoeHierarchical2.spkHistPrior3(params, results, spkHistPrior_setup);
    GMLMstructure.scaleParams          = @(params) DMC.priors.HorseshoeHierarchical2.scaleWB3(params, localCouplingGroup);
    GMLMstructure.scaleDerivatives     = @(results,params,posterior,A) DMC.priors.HorseshoeHierarchical2.scaleDWB3(results, params, localCouplingGroup, A);
    GMLMstructure.prior.dim_H          = spkHistPrior_setup.NH;
end

prior_setup = struct("V", cell(numel(GMLMstructure.Groups),1), "T", [], "tau", [], "phi", [], "lambda", [], "c", [], "log_constant_scale", []);
for jj = 1:numel(GMLMstructure.Groups)

    prior_setup(jj).log_constant_scale = log(modelSetup.constant_scale) + GMLMstructure.Groups(jj).extra_log_constant_scale;
    prior_setup(jj).include_rank_constant = true; % scale the prior over tau by the rank of the component sqrt(1/R) (sqrt because it's multiplied in twice - T{2} and V)
    prior_setup(jj).tau.log_scale = 0;
    prior_setup(jj).tau.dfs = 1;
    prior_setup(jj).phi.log_scale = 0;
    prior_setup(jj).phi.dfs = 1;
    prior_setup(jj).lambda.dfs = 1;
    prior_setup(jj).c.a = 2;  % t_{nu}(0, s^2)  => a=nu/2, b = nu*s^2/2  , nu = 2*a,  s = sqrt(2*b/nu) = sqrt(b/a), Piironen & Vehtari: a = 2, b = 8
    prior_setup(jj).c.b = 4;
    
    prior_setup(jj).V.on = true;
    prior_setup(jj).V.grps = [];
    prior_setup(jj).V.lambda_log_scale = GMLMstructure.Groups(jj).loadings_log_constant_scales;
    prior_setup(jj).V.mu = 0;
    prior_setup(jj).V.dim_T = GMLMstructure.dim_P;

    prior_setup(jj).T(1).on = false;
    prior_setup(jj).T(1).grps = [];
    prior_setup(jj).T(1).dim_T = GMLMstructure.Groups(jj).dim_T(1);
    prior_setup(jj).T(1).mu = 0;

    if(~GMLMstructure.Groups(jj).is_coupling)
        basis = lower(GMLMstructure.Groups(jj).name);
        prior_setup(jj).T(1).U_sigma_0         = prior_normalization.(basis).U_sigma_0;
        prior_setup(jj).T(1).U_sigma           = prior_normalization.(basis).U_sigma ;
        prior_setup(jj).T(1).U_transform       = prior_normalization.(basis).U_transform;
        prior_setup(jj).T(1).lambda_log_scale  = log(prior_normalization.(basis).sig);
        prior_setup(jj).T(1).U_sigma_inv       = inv(prior_setup(jj).T(1).U_sigma_0);
        prior_setup(jj).T(1).U_sigma_chol_inv  = chol(prior_setup(jj).T(1).U_sigma_inv)';
    else
        prior_setup(jj).T(1).U_sigma_0   = prior_normalization.spkHist.U_sigma_0;
        prior_setup(jj).T(1).U_sigma     = prior_normalization.spkHist.U_sigma ;
        prior_setup(jj).T(1).U_transform = prior_normalization.spkHist.U_transform;
        prior_setup(jj).T(1).lambda_log_scale = 0;
        prior_setup(jj).T(1).U_sigma_inv      = inv(prior_setup(jj).T(1).U_sigma_0);
        prior_setup(jj).T(1).U_sigma_chol_inv = chol(prior_setup(jj).T(1).U_sigma_inv)';
    end
    

    if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
        prior_setup(jj).T(2).on = true;
        prior_setup(jj).T(2).lambda_log_scale = GMLMstructure.Groups(jj).grp_log_constant_scales(:);
        prior_setup(jj).T(2).mu = 0;
        prior_setup(jj).T(2).dim_T     = GMLMstructure.Groups(jj).dim_T(2);
        prior_setup(jj).T(2).grps      = GMLMstructure.Groups(jj).grp_coef_idxs;
        prior_setup(jj).T(2).grp_names = GMLMstructure.Groups(jj).grp_coef_names;
    end

    if(numel(GMLMstructure.Groups(jj).dim_T) > 2)
        error("Prior not setup for this configuration");
    end

    if(~isfield(modelSetup, "flatPrior") || modelSetup.flatPrior)
        GMLMstructure.Groups(jj).prior.dim_H          = @(rank) DMC.priors.Horseshoe3.getDimH(rank, prior_setup(jj));
        GMLMstructure.Groups(jj).prior.log_prior_func = @(params, results, groupNum, priorOnly) DMC.priors.Horseshoe3.groupedRankHorseshoe(params, results, groupNum, prior_setup(jj), priorOnly);
        GMLMstructure.Groups(jj).prior.generator      = @(rank) DMC.priors.Horseshoe3.generateH(rank, prior_setup(jj), false, false);
    
        GMLMstructure.Groups(jj).scaleParams      = @(params) DMC.priors.Horseshoe3.scaleParams(params, prior_setup(jj));
        GMLMstructure.Groups(jj).scaleDerivatives = @(results,params,posterior,A) DMC.priors.Horseshoe3.scaleDerivatives(results, params, prior_setup(jj), posterior, A);
    else
        GMLMstructure.Groups(jj).prior.dim_H          = @(rank) DMC.priors.HorseshoeHierarchical.getDimH(rank, prior_setup(jj));
        GMLMstructure.Groups(jj).prior.log_prior_func = @(params, results, groupNum, priorOnly) DMC.priors.HorseshoeHierarchical.groupedRankHorseshoe(params, results, groupNum, prior_setup(jj));
        GMLMstructure.Groups(jj).prior.generator      = @(rank) DMC.priors.HorseshoeHierarchical.generateH(rank, prior_setup(jj), false, false);
    end
end
    

%% setup rescaling MH step
%parameters for an MH step to quickly traverse the scaler part of each component of tensor
if(modelSetup.useGibbsStep)
    MH_scaleSettings.sig =  0.2;
    MH_scaleSettings.N   = 10; % I used to sample 10 because I could, but 5 ought to be more than enough
    MH_scaleSettings.sample_every = 1;
    
    for jj = 1:numel(GMLMstructure.Groups)
        GMLMstructure.Groups(jj).gibbs_step.dim_H = 0;
        GMLMstructure.Groups(jj).gibbs_step.sample_func = @(gmlm, params, optStruct, sampleNum, groupNum, opts, results) DMC.priors.Horseshoe2.scalingMHStep(gmlm, params, optStruct, sampleNum, groupNum,  MH_scaleSettings, prior_setup(jj));
    end
else
    for jj = 1:numel(GMLMstructure.Groups)
        GMLMstructure.Groups(jj).gibbs_step.dim_H = 0;
        GMLMstructure.Groups(jj).gibbs_step.sample_func = [];
    end
end


end

