%function [LL_s, LL2_s, FR_s, dim_N_ranges, Ts_name, filterLength] = setupModelForEvaluation2(modelSetup_0, sample_idx)
%load('/home/latimerk/gitCode/GMLM_RCT/Results/GMLM/Coupling_local/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC32_LC_run1.mat_part.mat', 'modelInfo')
RCT.utils.setupGMLMpaths()
filterLength = 5;
pruneRate = 10;
if(~exist("modelInfo", "var") && exist("modelSetup", "var"))
    modelSetup_0 = modelSetup;
    sample_idxs = (8e3+1):pruneRate:(48e3);
else
    modelSetup_0 = modelInfo;
    sample_idxs = (HMC_settings.stepSize.scaleRanges(1)):pruneRate:(sample_idx-1);
end
if(isempty(sample_idxs))
    error("No samples");
end

folders = RCT.dataHandlers.getFolders();

locations = [modelSetup_0.location modelSetup_0.interAreaCoupling_locations];
downsample = round(modelSetup_0.delta_t_sec/1e-3);

[fname, fname_base, fname_samples, fname_samples_dat] = RCT.modelFitting.getModelFitFileName(modelSetup_0.subject, modelSetup_0.session, modelSetup_0);

fname_samples_partial1 = sprintf("%s_part1.mat", fname_samples);
fname_samples_partial2 = sprintf("%s_part2.mat", fname_samples);


if(exist(fname_samples, "file") || exist(fname_samples_partial1, "file") || exist(fname_samples_partial1, "file"))
    if(~exist(fname_samples, "file"))
        si_1.sample_idx = -1;
        si_2.sample_idx = -1;
        if(exist(fname_samples_partial1, "file"))
            si_1 = load(fname_samples_partial1, "sample_idx");
        end
        if(exist(fname_samples_partial2, "file"))
            si_2 = load(fname_samples_partial2, "sample_idx");
        end
        if(si_1.sample_idx > si_2.sample_idx)
            fname_samples = fname_samples_partial1;
        else
            fname_samples = fname_samples_partial2;
        end
    end

    %% load data
    data = RCT.dataHandlers.loadData(modelSetup_0.subject, modelSetup_0.session, locations, [], downsample);
    modelSetup = modelSetup_0;

    if(trainSet)
        qs = ismember(1:4, modelSetup.fitQuarter);
    else
        qs = ~ismember(1:4, modelSetup.fitQuarter);
    end
    if(all(~qs))
        modelSetup.fitQuarter = 1:4;
    else
        modelSetup.fitQuarter = find(qs);
    end
    modelSetup.includeTrialwiseLatent = false;
    modelSetup.interAreaCoupling_shuffle = false;
    
    %% get baseline rates
    if(~modelSetup.BASELINE_model && modelSetup.baselineCorrectionType ~= "none")
        baselineRates = RCT.modelBuilder.loadBaselineRates(modelSetup_0);
    else
        baselineRates = [];
    end

    %%
    [GMLMstructure, trials, ~, ~, trialsUsed] = RCT.modelBuilder.constructGMLMdata14(data, modelSetup, "baselineRates", baselineRates);
    gmlm = kgmlm.GMLM(GMLMstructure, trials, modelSetup.delta_t_sec, modelSetup.llType);

    %% get samples
    samples = RCT.resultsHandlers.loadSample(fname_samples, sample_idxs);
    [LL_s, LL2_s, FR_s, dim_N_ranges, Z] = RCT.resultsHandlers.evaluateLogLikelihood3(gmlm, samples, modelSetup_0, filterLength);

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