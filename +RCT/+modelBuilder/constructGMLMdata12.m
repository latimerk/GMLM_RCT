function [GMLMstructure, trials, spkHistPrior_setup, prior_setup, trialsUsed, latentBasis, couplingShuffle, normalizationSettings, timeWindows] = constructGMLMdata12(data, modelSetup, varargin)


p = inputParser;
p.CaseSensitive = false;

addRequired(p, "modelSetup"   ,    @(A) isstruct(A) );

default_noise_time = 0;
if(ismember("fixation_time", [modelSetup.stimConfig(:).event_name]))
    default_noise_time = 500;
end

addParameter(p, "t_pre_noise_bins"   ,  default_noise_time  / data.bin_size_ms,    @isnumeric);
addParameter(p, "t_post_saccade_bins",  60 / data.bin_size_ms,    @isnumeric);

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

% trQtrs = ceil(linspace(1, numel(trialsUsed)+1, 5));
% qtr1 = trialsUsed(trQtrs(1):(trQtrs(2)-1));
% qtr2 = trialsUsed(trQtrs(2):(trQtrs(3)-1));
% qtr3 = trialsUsed(trQtrs(3):(trQtrs(4)-1));
% qtr4 = trialsUsed(trQtrs(4):(trQtrs(5)-1));

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


%% Setup GMLM groups
GMLMstructure.dim_B = size(modelSetup.bases.spkHist.B, 2); % specify size of the linear term
GMLMstructure.dim_P = size(data.trials(1).Y.(modelSetup.location), 2);
GMLMstructure.Groups = struct("X_shared", [], "dim_R_max", [], "dim_A", [], "name", [], "dim_names", [], "dim_T", [], "dim_F", [], "factor_idx", [], "is_coupling", []); % tensor coefficient groups

dim_M = numel(data.trials);

%  stimulus setup (DENSE / local regressors)

events_stimulus = ["noise_on";
                    "stim_on";
                    "targets_on";
                    "fix_off";
                    "fixation_time"];
events_stimulus = events_stimulus(ismember(events_stimulus, [modelSetup.stimConfig(:).event_name]));




events_response = "saccade_end";



[~,event_idx] = ismember(events_stimulus, [modelSetup.stimConfig(:).event_name]);
[~,event_idx_r] = ismember(events_response, [modelSetup.responseConfig(:).event_name]);
BB = cell2mat({modelSetup.stimConfig(event_idx).regressors}');
cond_init = cumsum([0;arrayfun(@(a) size(a.regressors,1), modelSetup.stimConfig(event_idx) )]);

if(~modelSetup.includeTrialwiseLatent)
    jj = 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Stimuli";
    GMLMstructure.Groups(jj).dim_names = ["timing", "xstim"]; %dimensions of the tensor
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.stimulus; % max allocated space for rank
    GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.stimulus.B, 2) size(modelSetup.stimConfig(1).regressors,2)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    
    dense_stimulus = false;
    if(~dense_stimulus)
        GMLMstructure.Groups(jj).dim_A = numel(events_stimulus);    
        GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.stimulus.B, BB};
        GMLMstructure.Groups(jj).dim_F = (GMLMstructure.Groups(jj).dim_T);
        GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
    else
        GMLMstructure.Groups(jj).dim_A = 1;     
        GMLMstructure.Groups(jj).X_shared = {[]};
        GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.stimulus.B, 2) size(BB,2)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
        GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
        GMLMstructure.Groups(jj).factor_idx = [1 1]; %factor setup
    end
    
    % Response setup
    jj = jj + 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Response";
    GMLMstructure.Groups(jj).dim_A = 1;     
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.response; % max allocated space for rank
    if(size(modelSetup.responseConfig.regressors,2) > 1) % if a tensor is necessary
        dense_response = false;
    
        GMLMstructure.Groups(jj).dim_names = ["timing", "xchoice"]; %dimensions of the tensor
        GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.response.B, 2) size(modelSetup.responseConfig(1).regressors,2)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    
        if(~dense_response)
            GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.response.B, modelSetup.responseConfig(1).regressors};
            GMLMstructure.Groups(jj).dim_F = GMLMstructure.Groups(jj).dim_T;
            GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
        else
            GMLMstructure.Groups(jj).X_shared = {[]};
            GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
            GMLMstructure.Groups(jj).factor_idx = [1 1]; %factor setup
        end
    
    else %low rank matrix instead
    
        GMLMstructure.Groups(jj).dim_names = "timing"; %dimensions of the tensor
        GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.response.B };
        GMLMstructure.Groups(jj).dim_T = size(modelSetup.bases.response.B, 2) ; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
        GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
        GMLMstructure.Groups(jj).factor_idx = 1; %factor setup
    
        dense_response = false;
    end
else
    latentScaling = 1; %modelSetup.includeTrialwiseLatent_scaleFunc(dim_M)

    if(modelSetup.includeTrialwiseLatent_allStimEvents)
        stim_latent_groups = arrayfun(@(aa) size(aa.regressors,1), modelSetup.stimConfig(event_idx));
        stim_latent_groups_names = [modelSetup.stimConfig(event_idx).event_name];
    else
        stim_latent_groups = sum(arrayfun(@(aa) size(aa.regressors,1), modelSetup.stimConfig(event_idx)));
        stim_latent_groups_names = "stim";
    end

    jj = 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Stimuli";
    GMLMstructure.Groups(jj).dim_names = ["timing", "xstim_trial"]; %dimensions of the tensor
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.stimulus; % max allocated space for rank
    GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.stimulus.B, 2) size(BB,2)+dim_M*numel(stim_latent_groups)]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct

    MM = zeros(sum(stim_latent_groups) * dim_M, dim_M * numel(stim_latent_groups));
    for ii = 1:numel(stim_latent_groups)
        bb = zeros(sum(stim_latent_groups),1);
        bb((cond_init(ii)+1):(cond_init(ii+1))) = 1;

        MM(:, (1:dim_M) + (ii-1)*dim_M) = kron(latentScaling*eye(dim_M), bb);
    end
    
    dense_stimulus = false;
    if(~dense_stimulus)
        GMLMstructure.Groups(jj).dim_A = numel(events_stimulus);    
        GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.stimulus.B, [repmat(BB, dim_M, 1) MM]};
        GMLMstructure.Groups(jj).dim_F = (GMLMstructure.Groups(jj).dim_T);
        GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
    else
            error("not set up!");
        GMLMstructure.Groups(jj).dim_A = 1;     
        GMLMstructure.Groups(jj).X_shared = {[], eye(dim_M)};
        GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.stimulus.B, 2) size(modelSetup.stimConfig(1).regressors,2) dim_M]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
        GMLMstructure.Groups(jj).dim_F = [prod(GMLMstructure.Groups(jj).dim_T(1:2)) dim_M];
        GMLMstructure.Groups(jj).factor_idx = [1 1 2]; %factor setup
    end
    
    % Response setup
    jj = jj + 1;
    GMLMstructure.Groups(jj).is_coupling = false;
    GMLMstructure.Groups(jj).name = "Response";
    GMLMstructure.Groups(jj).dim_A = 1;     
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.response; % max allocated space for rank
    if(size(modelSetup.responseConfig.regressors,2) > 1) % if a tensor is necessary
        dense_response = false;
    
        GMLMstructure.Groups(jj).dim_names = ["timing", "xchoice_trial"]; %dimensions of the tensor
        BB_r = modelSetup.responseConfig(1).regressors;
        GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.response.B, 2) size(BB_r,2)+dim_M]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    
        if(~dense_response)
            GMLMstructure.Groups(jj).X_shared = {modelSetup.bases.response.B, [repmat(BB_r, dim_M, 1) kron(latentScaling * eye(dim_M), ones(size(BB_r,1),1))]};
            GMLMstructure.Groups(jj).dim_F = GMLMstructure.Groups(jj).dim_T;
            GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
        else
            error("not set up!");
            GMLMstructure.Groups(jj).X_shared = {[], eye(dim_M)};
            GMLMstructure.Groups(jj).dim_F = [prod(GMLMstructure.Groups(jj).dim_T(1:2)) dim_M];
            GMLMstructure.Groups(jj).factor_idx = [1 1 2]; %factor setup
        end
    
    else %low rank matrix instead
        error("not set up!");
    
        GMLMstructure.Groups(jj).dim_names = ["timing", "trial"]; %dimensions of the tensor
        GMLMstructure.Groups(jj).X_shared = {size(modelSetup.bases.response.B, 2), eye(dim_M) };
        GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.response.B, 2), dim_M] ; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
        GMLMstructure.Groups(jj).dim_F = GMLMstructure.Groups(jj).dim_T;
        GMLMstructure.Groups(jj).factor_idx = [1 2]; %factor setup
    
        dense_response = false;
    end
end



% Coupling
coupling_locations = [];
if(modelSetup.includeCouplingLocal)
    % local coupling
    jj = jj + 1;
    couplingShuffle.(modelSetup.location) = 1:dim_M;
    couplingShuffle.event.(modelSetup.location) = "noise_on";

    coupling_locations = cat(2, coupling_locations, modelSetup.location);
    GMLMstructure.Groups(jj).is_coupling = true;
    GMLMstructure.Groups(jj).name = "Coupling_Local";
    GMLMstructure.Groups(jj).dim_names = ["kernel", "neuron_weights"]; %dimensions of the tensor
    GMLMstructure.Groups(jj).dim_A = 1;     
    GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.coupling_local; % max allocated space for rank
    GMLMstructure.Groups(jj).X_shared = {[]};
    GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.spkCoupling.B, 2) GMLMstructure.dim_P]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
    GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
    GMLMstructure.Groups(jj).factor_idx = [1 1]; %factor setup

    NF_c = 1;

    if(modelSetup.includeCouplingLocal_trialDim)
        GMLMstructure.Groups(jj).dim_names  = cat(2, GMLMstructure.Groups(jj).dim_names, "trial_weights");
        GMLMstructure.Groups(jj).X_shared   = cat(2, GMLMstructure.Groups(jj).X_shared, eye(dim_M));
        GMLMstructure.Groups(jj).dim_T      = cat(2, GMLMstructure.Groups(jj).dim_T, dim_M);
        GMLMstructure.Groups(jj).dim_F      = cat(2, GMLMstructure.Groups(jj).dim_F, dim_M);
        GMLMstructure.Groups(jj).factor_idx = cat(2, GMLMstructure.Groups(jj).factor_idx, NF_c+1);
        NF_c = NF_c+1; %#ok<NASGU> 
    end
end
if(modelSetup.includeCouplingInter)
    % from other areas
    for ii = 1:numel(data.NeuronInfo)
        loc_c = data.NeuronInfo(ii).location;
        if(~strcmpi(loc_c, modelSetup.location))

            if(isfield(modelSetup, "interAreaCoupling_shuffle") && modelSetup.interAreaCoupling_shuffle)
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
                couplingShuffle.(loc_c) = 1:dim_M;
                couplingShuffle.event.(loc_c) = "noise_on";
            end

            jj = jj + 1;
            coupling_locations = cat(2, coupling_locations, loc_c);
            GMLMstructure.Groups(jj).is_coupling = true;
            dim_P_c = size(data.trials(1).Y.(loc_c),2);
            GMLMstructure.Groups(jj).name = sprintf("Coupling_%s", loc_c);
            GMLMstructure.Groups(jj).dim_names = ["kernel", "neuron_weights"]; %dimensions of the tensor
            GMLMstructure.Groups(jj).dim_A = 1;     
            GMLMstructure.Groups(jj).dim_R_max = modelSetup.Ranks.coupling_inter; % max allocated space for rank
            GMLMstructure.Groups(jj).X_shared = {[]};
            GMLMstructure.Groups(jj).dim_T = [size(modelSetup.bases.spkCoupling.B, 2) dim_P_c]; % I require specifying the dimensions of each part of the tensor to make sure everything is correct
            GMLMstructure.Groups(jj).dim_F = prod(GMLMstructure.Groups(jj).dim_T);
            GMLMstructure.Groups(jj).factor_idx = [1 1]; %factor setup

            NF_c = 1;
        
            if(modelSetup.includeCouplingInter_trialDim)
                GMLMstructure.Groups(jj).dim_names  = cat(2, GMLMstructure.Groups(jj).dim_names, "trial_weights");
                GMLMstructure.Groups(jj).X_shared   = cat(2, GMLMstructure.Groups(jj).X_shared, eye(dim_M));
                GMLMstructure.Groups(jj).dim_T      = cat(2, GMLMstructure.Groups(jj).dim_T, dim_M);
                GMLMstructure.Groups(jj).dim_F      = cat(2, GMLMstructure.Groups(jj).dim_F, dim_M);
                GMLMstructure.Groups(jj).factor_idx = cat(2, GMLMstructure.Groups(jj).factor_idx, NF_c+1);
                NF_c = NF_c+1; %#ok<NASGU> 
            end
        end
    end
end
N_tensor_groups = numel(GMLMstructure.Groups);

%% sets up trials
trials = struct("Y", cell(numel(data.trials),1), "X_lin", [], "neuron", [], "Groups", []);

paddedHspk     = [zeros(size(modelSetup.bases.spkHist.B)     + [1 0]); modelSetup.bases.spkHist.B];%zero padding is to make convolutions easier
paddedCoupling = [zeros(size(modelSetup.bases.spkCoupling.B) + [1 0]); modelSetup.bases.spkCoupling.B];

N_bins_response = 0;
N_bins_stim     = 0;
N_bins_trial    = nan(numel(data.trials),1);
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
    trials(tt).X_lin = zeros([TT size(modelSetup.bases.spkHist.B,2) GMLMstructure.dim_P]); 

    Y_t = data.trials(tt).Y.(modelSetup.location);

    for bb = 1:size(modelSetup.bases.spkHist.B,2)
        h_c = conv2(Y_t, paddedHspk(:, bb), "same");
        trials(tt).X_lin(:, bb, :) = reshape(h_c(tts_c, :), TT,1,[]);
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

        tt_c = couplingShuffle.(loc_c)(tt);
        dt_0_coupling = data.trials(tt_c).(couplingShuffle.event.(loc_c) ) - data.trials(tt).(couplingShuffle.event.(loc_c) );
        tts_c_coupling = tts_c + dt_0_coupling;
        if(tts_c_coupling(1) < 1 || tts_c_coupling(end) > size(data.trials(tt_c).Y.(loc_c),1))
            error("Trial window exceeds processed data window size!");
        end

        dim_P_c = size(data.trials(tt_c).Y.(loc_c),2);

        NF_c = numel(unique(GMLMstructure.Groups(jj).factor_idx));
        trials(tt).Groups(jj).X_local    = cell(1,NF_c); 
        trials(tt).Groups(jj).iX_shared  = cell(1,NF_c); 
        trials(tt).Groups(jj).X_local{1} = zeros([TT size(modelSetup.bases.spkCoupling.B,2) dim_P_c]);  

        Y_t = data.trials(tt_c).Y.(loc_c);
        

        for bb = 1:size(modelSetup.bases.spkCoupling.B,2)
            h_c = conv2(Y_t, paddedCoupling(:, bb), "same");
            trials(tt).Groups(jj).X_local{1}(:, bb, :) = reshape(h_c(tts_c_coupling, :), TT, 1, []);
        end  

        ff = strcmpi(GMLMstructure.Groups(jj).dim_names, "trial_weights");
        if(~all(~ff))
            ff = GMLMstructure.Groups(jj).factor_idx(ff);
            trials(tt).Groups(jj).iX_shared{ff} = ones(TT,1) * tt;
        end 
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
                  target_config;%fix_off
                  1]; %fixation time
    
    jj = 1;
    X = zeros([TT GMLMstructure.Groups(jj).dim_F]);  
    trials(tt).Groups(jj).iX_shared = {[]}; 
    if(dense_stimulus)
        X   = zeros([TT GMLMstructure.Groups(jj).dim_F]); 
    else
        iX  = zeros(TT, numel(events_stimulus)); 
        iX2 = zeros(TT, numel(events_stimulus)); 
    end

    for aa = 1:numel(events_stimulus)
        timing = data.trials(tt).(events_stimulus(aa)) - t_0 + 1;
        if(isnan(timing))
            continue;
        end

        tts_event = timing + modelSetup.bases.stimulus.tts;
        vv = tts_event > 0 & tts_event <= TT;
        if(all(~vv))
            continue;
        end

        if(dense_stimulus)
            B = kron(modelSetup.stimConfig(event_idx(aa)).regressors(conditions(aa),:), modelSetup.bases.stimulus.B); %#ok<*UNRCH> 
            X(tts_event(vv),:) = X(tts_event(vv),:) + B(vv,:);
        else
            iX( tts_event(vv), aa) = modelSetup.bases.stimulus.tts_idx(vv); 
            iX2(tts_event(vv), aa) = conditions(aa) + cond_init(event_idx(aa));
        end
    end
    if(~modelSetup.includeTrialwiseLatent)  
        if(dense_stimulus)
            trials(tt).Groups(jj).iX_shared = {[]};
            trials(tt).Groups(jj).X_local = {X};
        else
            trials(tt).Groups(jj).iX_shared = {iX, iX2};
            trials(tt).Groups(jj).X_local = {[], []};
        end
    else  
        if(dense_stimulus)
            trials(tt).Groups(jj).iX_shared = {[], ones(TT,GMLMstructure.Groups(jj).dim_A)*tt};
            trials(tt).Groups(jj).X_local = {X, []};
        else
            pp = ones(size(iX2)) * size(BB,1)*(tt-1);
            pp(iX2 <= 0) = 0;
            trials(tt).Groups(jj).iX_shared = {iX, iX2 + pp};
            trials(tt).Groups(jj).X_local = {[], []};
        end
    end

    %% setup response
    conditions = target_choice_color + (target_choice_location-1)*2;

    jj = 2;

    if(dense_response)
        X = zeros([TT GMLMstructure.Groups(jj).dim_F]);  
        iX = []; 
        iX2 = [];
    else
        X = [];
        iX  = zeros(TT, numel(events_response)); 
        if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
            iX2 = zeros(TT, numel(events_response)); 
        end
    end
    X2 = [];
    for aa = 1:numel(events_response)
        timing = data.trials(tt).(events_response(aa)) - t_0 + 1;
        if(isnan(timing))
            continue;
        end

        tts_event = timing + modelSetup.bases.response.tts;
        vv = tts_event > 0 & tts_event <= TT;

        if(all(~vv))
            continue;
        end


        if(dense_response)
            B = kron(modelSetup.responseConfig(event_idx_r(aa)).regressors(conditions(aa),:), modelSetup.bases.response.B);
            X(tts_event(vv),:) = X(tts_event(vv),:) + B(vv,:);
        else
            iX(tts_event(vv), aa) = modelSetup.bases.response.tts_idx(vv); %#ok<AGROW>
            if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
                iX2(tts_event(vv), aa) = conditions(aa);
            end
        end

        if(aa == 1)
            % get num bins and mean spike count contributing to response window
            N_bins_response = N_bins_response + numel(tts_event(vv));
            Spk_count_response = Spk_count_response + sum(trials(tt).Y(tts_event(vv),:),1);
        end
    end
    if(~modelSetup.includeTrialwiseLatent)  
        if(numel(GMLMstructure.Groups(jj).dim_F) > 1)
            trials(tt).Groups(jj).iX_shared = {iX, iX2};
            trials(tt).Groups(jj).X_local = {X, X2};
        else
            trials(tt).Groups(jj).iX_shared = {iX};
            trials(tt).Groups(jj).X_local = {X};
        end
    else
        if(numel(GMLMstructure.Groups(jj).dim_F) > 1)
            pp = ones(size(iX2)) * size(BB_r,1)*(tt-1);
            pp(iX2 <= 0) = 0;
            trials(tt).Groups(jj).iX_shared = {iX, iX2 + pp};
            trials(tt).Groups(jj).X_local = {X, X2};
        else
            trials(tt).Groups(jj).iX_shared = {iX};
            trials(tt).Groups(jj).X_local = {X};
        end

    end
end





%%
fprintf("Done.\n");




%% setup Horseshoe priors to encourage sparsity
% f_0   = repmat(1/4, numel(GMLMstructure.Groups));
% mu_stim = mean(Spk_count_stim    )/N_bins_stim;
% mu_resp = mean(Spk_count_response)/N_bins_response;

% if(modelSetup.includeTrialwiseLatent)
%     N_bins_trial_c = mean(N_bins_trial);
%     mu_trial = mean(Spk_count_trial,"all")/N_bins_trial_c;
%     M = [N_bins_stim; N_bins_response; N_bins_trial_c; repmat(N_bins_stim, N_tensor_groups-3, 1)] * GMLMstructure.dim_P;
%     mu_tilde = [mu_stim; mu_resp; mu_trial; repmat(mu_stim, N_tensor_groups-3, 1)] ;
% else
%     M = [N_bins_stim; N_bins_response; repmat(N_bins_stim, N_tensor_groups-2, 1)] * GMLMstructure.dim_P;
%     mu_tilde = [mu_stim; mu_resp; repmat(mu_stim, N_tensor_groups-2, 1)] ;
% end

coef_lambda_scales = {1; 1; 1; 1; 1; 1};
for jj = 3:numel(GMLMstructure.Groups)
    coef_lambda_scales{jj} = 1;%/sqrt(GMLMstructure.Groups(jj).dim_T(2));
end

if(~modelSetup.includeTrialwiseLatent)
    coef_names = {[modelSetup.stimConfig(:).coeff_names]; [modelSetup.responseConfig(:).coeff_names]};
    coef_idxs = {[modelSetup.stimConfig(:).coeff_groups]; [modelSetup.responseConfig(:).coeff_groups]};
else
    cn_stim = [modelSetup.stimConfig(:).coeff_names];
    ci_stim = [modelSetup.stimConfig(:).coeff_groups];
    for ii = 1:numel(stim_latent_groups)
        cn_stim = [cn_stim sprintf("trial_%s", stim_latent_groups_names(ii))]; %#ok<AGROW> 
        ci_stim = [ci_stim (size(BB,2)+ (ii-1)*dim_M + (1:dim_M) )]; %#ok<AGROW> 
    end

    coef_names = {cn_stim [modelSetup.responseConfig(:).coeff_names  "trial_response"]};
    coef_idxs = {ci_stim; [modelSetup.responseConfig(:).coeff_groups (size(BB_r,2) + (1:dim_M))]};

end

%% removes unused latents
if(modelSetup.reprojectLatents)
    rng_0 = rng();
    rng(20220501);
end
if(modelSetup.includeTrialwiseLatent)
    for jj = 1:2
        if(jj == 1) 
            NM = size(BB,2);
        else
            NM = size(BB_r,2);
        end
        NA = GMLMstructure.Groups(jj).dim_A;
        cs = true(NM + dim_M*NA,1);


        %% finds any latents for events that didn't occur
        if(modelSetup.removeUnusedLatents)
            for mm = 1:dim_M
                unused_idx = find(all(trials(mm).Groups(jj).iX_shared{2} == 0,1));
                if(~isempty(unused_idx))
                    unused_idx = (unused_idx-1)*dim_M + mm + NM;
                    cs(unused_idx) = false;
                end
            end
    
            if(~all(cs))
                %% remove columns from shared idx
                GMLMstructure.Groups(jj).X_shared{2} = GMLMstructure.Groups(jj).X_shared{2}(:, cs);
                %% remove columns from coef_idxs
                aa = find(cs);
                for gg = 1:numel(coef_idxs{jj})
                    [~,cc] = ismember(coef_idxs{jj}{gg}, aa);
                    coef_idxs{jj}{gg} = cc(cc ~= 0);
                end
    
                GMLMstructure.Groups(jj).dim_T(2) = size(GMLMstructure.Groups(jj).X_shared{2},2);
                GMLMstructure.Groups(jj).dim_F(2) = size(GMLMstructure.Groups(jj).X_shared{2},2);
    
            end
        end
    
        coef_lambda_scales{jj} = ones(numel(coef_idxs{jj}), 1);
        for gg = 1:numel(coef_idxs{jj})
            if(startsWith(coef_names{jj}(gg), "trial"))

                ls = modelSetup.includeTrialwiseLatent_scaleFunc(dim_M);
%                 ls = modelSetup.includeTrialwiseLatent_scaleFunc(numel(coef_idxs{jj}{gg}));
                coef_lambda_scales{jj}(gg) = ls;
            end
        end

        %% re-organize basis: use orthonormal projection
        lb = cell(GMLMstructure.Groups(jj).dim_A,1);
        if(modelSetup.reprojectLatents)
            gg_ctr = 1;
            for gg = 1:numel(coef_idxs{jj})
                if(startsWith(coef_names{jj}(gg), "trial"))
                    ii = coef_idxs{jj}{gg} ;
                    X = GMLMstructure.Groups(jj).X_shared{2}(:, ii);
    
                    rng_0 = rng();
                    rng(gg + 1000*jj);
                    B = orth(randn(size(X,2), size(X,2)));
                    rng(rng_0);
                    Z = X*B;
    
                    GMLMstructure.Groups(jj).X_shared{2}(:, ii) = Z;
    
                    lb{gg_ctr} = B;
                    gg_ctr = gg_ctr + 1;
                end
            end
        else
            gg_ctr = 1;
            for gg = 1:numel(coef_idxs{jj})
                if(startsWith(coef_names{jj}(gg), "trial"))
                    ii = coef_idxs{jj}{gg};
                    X = GMLMstructure.Groups(jj).X_shared{2}(:, ii);
                    lb{gg_ctr} = eye(size(X,2));
                    gg_ctr = gg_ctr + 1;
                end
            end
        end
        
        if(jj == 1)
            latentBasis.stimulus = lb;
        else
            latentBasis.response = lb;
        end

    end
else
    latentBasis.stimulus = [];
    latentBasis.response = [];
end
if(modelSetup.reprojectLatents)
    rng(rng_0);
end
%%
    


if(~isfield(modelSetup, "useRidge") || ~modelSetup.useRidge)
    if(modelSetup.includeCouplingLocal)
        couplingGroup = 3;
    else
        couplingGroup = [];
    end

    spkHistPrior_setup.hyperprior.w_sig_nu = 1;
    spkHistPrior_setup.hyperprior.b_sig_nu = 1;
    spkHistPrior_setup.hyperprior.log_w_sig_scale = log(2);
    spkHistPrior_setup.hyperprior.log_b_sig_scale = log(2);
    spkHistPrior_setup.hyperprior.log_w_mu_sig = log(4);
    spkHistPrior_setup.hyperprior.log_b_mu_sig = log(4);

    if(GMLMstructure.dim_B > 0)
        spkHistPrior_setup.NH              = 2 + 1 + GMLMstructure.dim_B; %number of hyperparameters
    else
        spkHistPrior_setup.NH              = 2; %number of hyperparameters
    end

    if(modelSetup.kernelTimescales.spkHist > 0 && GMLMstructure.dim_B > 0)
        B = modelSetup.bases.spkHist.B;
        L = modelSetup.kernelTimescales.spkHist;
        tts = modelSetup.bases.spkHist.tts_0;
        stdParams = modelSetup.bases.spkHist.std;
        D = kernelFunction(tts, L, stdParams);
        S = B'*(D*B);

        spkHistPrior_setup.U_sigma_0 = S;
        spkHistPrior_setup.U_sigma = eye(size(S));
        spkHistPrior_setup.U_transform = chol(S);
    end
    
    GMLMstructure.prior.log_prior_func = @(params, results, priorOnly) DMC.priors.Horseshoe3.spkHistPrior3(params, results, spkHistPrior_setup, priorOnly);
    GMLMstructure.scaleParams      = @(params) DMC.priors.Horseshoe3.scaleWB3(params, spkHistPrior_setup, couplingGroup);
    GMLMstructure.scaleDerivatives = @(results,params,posterior,A) DMC.priors.Horseshoe3.scaleDWB3(results, params, spkHistPrior_setup, posterior, couplingGroup, A);

    
    GMLMstructure.prior.dim_H          = spkHistPrior_setup.NH;

    prior_setup = struct("V", cell(numel(GMLMstructure.Groups),1), "T", [], "tau", [], "phi", [], "lambda", [], "c", [], "log_constant_scale", []);
    for jj = 1:numel(GMLMstructure.Groups)
    
        prior_setup(jj).log_constant_scale = log(modelSetup.constant_scale);
        prior_setup(jj).include_rank_constant = true;
        prior_setup(jj).tau.log_scale = 0;
        prior_setup(jj).tau.dfs = 1;
        prior_setup(jj).phi.log_scale = 0;
        prior_setup(jj).phi.dfs = 1;
        prior_setup(jj).lambda.dfs = 1;
        prior_setup(jj).c.a = 2;  % t_{nu}(0, s^2)  => a=nu/2, b = nu*s^2/2  , nu = 2*a,  s = sqrt(2*b/nu) = sqrt(b/a), Piironen & Vehtari: a = 2, b = 8
        prior_setup(jj).c.b = 4;
        
        prior_setup(jj).V.on = true;
        prior_setup(jj).V.grps = [];
        prior_setup(jj).V.lambda_log_scale = 0;
        prior_setup(jj).V.mu = 0;
        prior_setup(jj).V.dim_T = GMLMstructure.dim_P;
    
    
        prior_setup(jj).T(1).on = false;
        prior_setup(jj).T(1).grps = [];
        prior_setup(jj).T(1).lambda_log_scale = 0;
        prior_setup(jj).T(1).dim_T = GMLMstructure.Groups(jj).dim_T(1);


        prior_setup(jj).T(1).mu = 0;
    
        if(jj > 2)
            B = modelSetup.bases.spkCoupling.B;
            L = modelSetup.kernelTimescales.spkCoupling;
            stdParams = modelSetup.bases.spkCoupling.std;
            tts = modelSetup.bases.spkCoupling.tts_0;
        elseif(jj == 2)
            B = modelSetup.bases.response.B;
            L = modelSetup.kernelTimescales.response;
            stdParams = modelSetup.bases.response.std;
            tts = modelSetup.bases.response.tts_0;
        else
            B = modelSetup.bases.stimulus.B;
            L = modelSetup.kernelTimescales.stimulus;
            stdParams = modelSetup.bases.stimulus.std;
            tts = modelSetup.bases.stimulus.tts_0;
        end
        
        if(L > 0)
            D = kernelFunction(tts, L, stdParams);
            S = B'*(D*B);
            prior_setup(jj).T(1).U_sigma_0 = S;
            prior_setup(jj).T(1).U_sigma = eye(size(S));
            prior_setup(jj).T(1).U_transform = chol(S);
        end

        if(numel(GMLMstructure.Groups(jj).dim_T) > 1)
            prior_setup(jj).T(2).on = true;
            prior_setup(jj).T(2).lambda_log_scale = log(coef_lambda_scales{jj});
            prior_setup(jj).T(2).mu = 0;
            prior_setup(jj).T(2).dim_T = GMLMstructure.Groups(jj).dim_T(2);
            if(jj <= numel(coef_idxs))
                prior_setup(jj).T(2).grps = coef_idxs{jj};
                prior_setup(jj).T(2).grp_names = coef_names{jj};
            else
                prior_setup(jj).T(2).grps = [];
                prior_setup(jj).T(2).grp_names = [];
            end
        end
    
    
        if(numel(GMLMstructure.Groups(jj).dim_T) > 2)
            error("Prior not setup for this configuration");
        end
    
        GMLMstructure.Groups(jj).prior.dim_H          = @(rank) DMC.priors.Horseshoe3.getDimH(rank, prior_setup(jj));
        GMLMstructure.Groups(jj).prior.log_prior_func = @(params, results, groupNum, priorOnly) DMC.priors.Horseshoe3.groupedRankHorseshoe(params, results, groupNum, prior_setup(jj), priorOnly);
        GMLMstructure.Groups(jj).prior.generator      = @(rank) DMC.priors.Horseshoe3.generateH(rank, prior_setup(jj), false, false);
    
        GMLMstructure.Groups(jj).scaleParams      = @(params) DMC.priors.Horseshoe3.scaleParams(params, prior_setup(jj));
        GMLMstructure.Groups(jj).scaleDerivatives = @(results,params,posterior,A) DMC.priors.Horseshoe3.scaleDerivatives(results, params, prior_setup(jj), posterior, A);
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
else
    warning("Using ridge prior!");
    spkHistPrior_setup.hyperprior.nu = 3;
    spkHistPrior_setup.NH            = 2; %number of hyperparameters
    GMLMstructure.prior.log_prior_func = @(params, results, priorOnly) DMC.priors.GMLMprior_spkHistConst(params, results, spkHistPrior_setup);
    GMLMstructure.prior.dim_H          = spkHistPrior_setup.NH;

    prior_setup = struct("V", cell(numel(GMLMstructure.Groups),1), "T", [], "tau", [], "phi", [], "lambda", [], "c", [], "log_constant_scale", []);
    for jj = 1:numel(GMLMstructure.Groups)
    
        prior_setup(jj).tau.log_scale = 0;
        prior_setup(jj).tau.dfs = 1;
        
        prior_setup(jj).V.on = true;
        prior_setup(jj).V.mu = 0;
    
        prior_setup(jj).T(1).mu = 0;
   
        prior_setup(jj).T(2).on = true;
        prior_setup(jj).T(2).mu = 0;
    
        GMLMstructure.Groups(jj).prior.dim_H          = @(rank) DMC.priors.Ridge.getDimH(rank, prior_setup(jj));
        GMLMstructure.Groups(jj).prior.log_prior_func = @(params, results, groupNum, priorOnly) DMC.priors.Ridge.Ridge(params, results, groupNum, prior_setup(jj), priorOnly);
        GMLMstructure.Groups(jj).prior.generator      = @(rank) DMC.priors.Ridge.generateH(rank, prior_setup(jj), false, false);
    
        GMLMstructure.Groups(jj).scaleParams      = @(params) DMC.priors.Ridge.scaleParams(params, prior_setup(jj));
        GMLMstructure.Groups(jj).scaleDerivatives = @(results,params,posterior,A) DMC.priors.Ridge.scaleDerivatives(results, params, prior_setup(jj), posterior, A);
    end

end





%% Scales any coupling and linear terms
if(~isfield(modelSetup, "normalizationSettings"))
    %% normalizes stim & response kernels
    fprintf("Normalizing stimulus & response regressors.\n");
    for jj = 1:2
        B = GMLMstructure.Groups(jj).X_shared{1};
        cS = prior_setup(jj).T(1).U_transform;
        C  = sqrt(trace(cov(B * cS))) * sqrt(2/pi);
        modelSetup.normalizationSettings.Groups(jj).std = C;
        modelSetup.normalizationSettings.Groups(jj).mu = 0;
    end

    if(GMLMstructure.dim_B > 0)
        fprintf("Estimating mean rates to subtract for spike history/coupling.\n");
        xx = cell2mat(arrayfun(@(aa)aa.X_lin, trials, 'UniformOutput', false));
        mu = mean(xx);
        modelSetup.normalizationSettings.spkHist.mu = mu;
        modelSetup.normalizationSettings.spkHist.std = zeros([1,1,size(xx,3)]);
        cS = spkHistPrior_setup.U_transform;
        for pp = 1:size(xx,3)
            modelSetup.normalizationSettings.spkHist.std(pp) = sqrt(trace(cov(xx(:,:,pp) * cS))) * sqrt(2/pi);
        end
    else
        modelSetup.normalizationSettings.spkHist.mu = [];
        modelSetup.normalizationSettings.spkHist.std = [];
    end

    for jj = 3:numel(GMLMstructure.Groups)
        cS = prior_setup(jj).T(1).U_transform;
        if(GMLMstructure.Groups(jj).name == "Coupling_Local")
            modelSetup.normalizationSettings.Groups(jj).mu = modelSetup.normalizationSettings.spkHist.mu;
            modelSetup.normalizationSettings.Groups(jj).std = modelSetup.normalizationSettings.spkHist.std;
        else
            xx = cell2mat(arrayfun(@(aa)aa.Groups(jj).X_local{1}, trials, 'UniformOutput', false));
            mu = mean(xx);
            modelSetup.normalizationSettings.Groups(jj).mu = mu;
            modelSetup.normalizationSettings.Groups(jj).std = zeros([1,1,size(xx,3)]);
            for pp = 1:size(xx,3)
                modelSetup.normalizationSettings.Groups(jj).std(pp) = sqrt(trace(cov(xx(:,:,pp) * cS)))* sqrt(2/pi);
            end
        end
    end
else
    fprintf("Using provided normalization settings.\n");
end


for jj = 1:2
    B = GMLMstructure.Groups(jj).X_shared{1};
    GMLMstructure.Groups(jj).X_shared{1} = (B - modelSetup.normalizationSettings.Groups(jj).mu) ./ modelSetup.normalizationSettings.Groups(jj).std;
end

for mm = 1:dim_M
    if(GMLMstructure.dim_B > 0)
        trials(mm).X_lin = (trials(mm).X_lin-modelSetup.normalizationSettings.spkHist.mu)./modelSetup.normalizationSettings.spkHist.std;
    end
    for jj = 3:numel(GMLMstructure.Groups)
        trials(mm).Groups(jj).X_local{1} = (trials(mm).Groups(jj).X_local{1} - modelSetup.normalizationSettings.Groups(jj).mu) ./ modelSetup.normalizationSettings.Groups(jj).std;
    end
end
normalizationSettings = modelSetup.normalizationSettings;

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