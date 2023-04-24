%function [] = setupModelForEvaluation(modelSetup_0, sample_idx)
%load('/home/latimerk/gitCode/GMLM_RCT/Results/GMLM/Coupling_local/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC32_LC_run1.mat_part.mat', 'modelInfo')
modelSetup_0 = modelInfo;
sample_idxs = (HMC_settings.nWarmup+1):200:(sample_idx-1);

folders = RCT.dataHandlers.getFolders();

locations = [modelSetup_0.location modelSetup_0.interAreaCoupling_locations];
downsample = round(modelSetup_0.delta_t_sec/1e-3);

[fname, fname_base, fname_samples, fname_samples_dat] = RCT.modelFitting.getModelFitFileName(modelSetup_0.subject, modelSetup_0.session, modelSetup_0);

fname_samples_partial = sprintf("%s_part.mat", fname_samples);


if(exist(fname_samples, "file") || exist(fname_samples_partial, "file"))
    if(~exist(fname_samples, "file"))
        fname_samples = fname_samples_partial;
    end

    %% load data
    data = RCT.dataHandlers.loadData(modelSetup_0.subject, modelSetup_0.session, locations, [], downsample);
    modelSetup = modelSetup_0;
    qs = ~ismember(1:4, modelSetup.fitQuarter);
    if(all(~qs))
        modelSetup.fitQuarter = 1:4;
    else
        modelSetup.fitQuarter = find(qs);
    end
    modelSetup.includeTrialwiseLatent = false;
    modelSetup.interAreaCoupling_shuffle = false;
    
    [GMLMstructure, trials, ~, trialsUsed] = RCT.modelBuilder.constructGMLMdata9(data, modelSetup);
    gmlm = kgmlm.GMLM(GMLMstructure, trials, modelSetup.delta_t_sec, modelSetup.llType);

    %% get samples
    samples = RCT.resultsHandlers.loadSample(fname_samples, sample_idxs);
    [LL, LL_prev, LL_curr, SpikeRate] = RCT.resultsHandlers.evaluateLogLikelihood(gmlm, samples, modelSetup_0);

    %% get trial timings
    Ts_name = ["noise_on", "stim_on", "targets_on", "fix_off", "saccade_end"];
    for mm = 1:numel(trialsUsed)
        for ff = 1:numel(Ts_name)
            if(isempty(data.trials(trialsUsed(mm)).(Ts_name(ff))))
                data.trials(trialsUsed(mm)).(Ts_name(ff)) = nan;
            end
        end
    end

    T_start = [data.trials(trialsUsed).noise_on] - 1;

    T_noise   = [data.trials(trialsUsed).noise_on];
    T_stim    = [data.trials(trialsUsed).stim_on];
    T_targets = [data.trials(trialsUsed).targets_on];
    T_fixoff  = [data.trials(trialsUsed).fix_off];
    T_sacc    = [data.trials(trialsUsed).saccade_end];

    Ts = [T_noise(:) T_stim(:) T_targets(:) T_fixoff(:) T_sacc(:)] - T_start(:);
else
    error("No results found!");
end