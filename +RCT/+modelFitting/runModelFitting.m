function [] = runModelFitting(data, GPUs,  varargin)

RCT.utils.setupGMLMpaths();

p = inputParser;
p.CaseSensitive = false;

addRequired(p, "GPUs",     @(A) isnumeric(A)  && ~isempty(A));
addOptional(p, "locations",    "", @(A) ischar(A) || isstring(A) || isempty(A));
p = RCT.modelFitting.setupDefaultModelFittingArgs(p);


parse(p, GPUs, varargin{:});
% then set/get all the inputs out of this structure
GPUs  = p.Results.GPUs;
locations  = p.Results.locations;
modelSetup.Ranks                  = p.Results.Ranks;
modelSetup.includeCouplingLocal   = p.Results.includeCouplingLocal;
modelSetup.includeCouplingInter   = p.Results.includeCouplingInter;
modelSetup.includeTrialwiseLatent = p.Results.includeTrialwiseLatent;
modelSetup.includeTrialwiseLatent_allStimEvents = p.Results.includeTrialwiseLatent_allStimEvents;
modelSetup.includeTrialwiseLatent_log_scaleFunc = p.Results.includeTrialwiseLatent_log_scaleFunc;
modelSetup.constant_scale = 0.1;

modelSetup.baselineCorrectionType = p.Results.baselineCorrectionType;

modelSetup.fitHalf                = p.Results.fitHalf;
modelSetup.fitQuarter             = p.Results.fitQuarter;
if(all(ismember(1:4, modelSetup.fitQuarter), "all"))
    modelSetup.fitQuarter = 0;
end

modelSetup.includeCouplingInter_trialDim = false; 
modelSetup.includeCouplingLocal_trialDim = false;

modelSetup.runNumber              = p.Results.runNumber;
modelSetup.DEBUG                  = p.Results.DEBUG;
modelSetup.doublePrecision        = p.Results.doublePrecision;
modelSetup.saveSinglePrecision    = p.Results.saveSinglePrecision;
modelSetup.useGibbsStep = p.Results.useGibbsStep;
modelSetup.BASELINE_model = p.Results.useBaseline;

modelSetup.reprojectLatents    = false;
modelSetup.removeUnusedLatents = false;

modelSetup.interAreaCoupling_shuffle = p.Results.shuffleCoupling;
modelSetup.NULL_Model = p.Results.useNull;


hmc_steps = p.Results.maxHMCsteps;
hmc_steps_0 = p.Results.maxHMCsteps_0;
HMCstep_L = p.Results.HMCstep_L;
HMC_delta = p.Results.HMC_delta;
HMC_delta_0 = p.Results.HMC_delta_0;
nSamples = p.Results.nSamples;
nWarmup = p.Results.nWarmup;


showPlots = p.Results.showPlots;
saveDuring = p.Results.saveDuring;
[~,IT]  = unique(p.Results.targets_verical_configs);
targets_verical_configs = p.Results.targets_verical_configs(sort(IT));
overwrite  = p.Results.overwrite;

validate_data = @(A) isstruct(A) && isfield(A, "NeuronInfo") && isfield(A,"TaskInfo") && isfield(A, "trials") && isstruct(A.trials) && isfield(A.trials,"Y");
if(~validate_data(data))
    error("Invalid data struct");
end

if(saveDuring && showPlots)
    error("Cannot currently plot progress while saving out samples.");
end

all_locations = [data.NeuronInfo(:).location];
if(isempty(locations))
    locations = all_locations;
end
modelSetup.includeCouplingInter = modelSetup.includeCouplingInter && numel(all_locations) > 1;


if(modelSetup.BASELINE_model && (modelSetup.NULL_Model || modelSetup.includeCouplingInter || modelSetup.includeCouplingLocal || modelSetup.includeTrialwiseLatent))
    error("Invalid options for baseline model!");
end


modelSetup.subject   = data.subject;
modelSetup.session   = data.session;
%% get model setup information
modelSetup.bases = RCT.modelBuilder.setupBasis(data.bin_size_ms);
modelSetup.kernelTimescales =  RCT.modelBuilder.getDefaultKernelTimescales();
if(isfield(modelSetup, "BASELINE_model") && modelSetup.BASELINE_model)
    [modelSetup.stimConfig, modelSetup.responseConfig, modelSetup.fixationConfig] = RCT.modelBuilder.getTaskBaselineModelSetup(data.TaskInfo);
    modelSetup.bases.spkHist.B = [];
elseif(isfield(modelSetup, "NULL_Model") && modelSetup.NULL_Model)
    [modelSetup.stimConfig, modelSetup.responseConfig] = RCT.modelBuilder.getTaskNullModelSetup(data.TaskInfo);
else
    [modelSetup.stimConfig, modelSetup.responseConfig] = RCT.modelBuilder.getTaskModelSetup(data.TaskInfo);
end



%%
for loc_idx = 1:numel(locations)
    modelSetup.location = locations(loc_idx);

    modelSetup.interAreaCoupling_locations = [];
    if(modelSetup.includeCouplingInter)
        modelSetup.interAreaCoupling_locations = all_locations(all_locations ~= modelSetup.location);
    end

    for targ_idx = 1:numel(targets_verical_configs)
        modelSetup.targets_verical             = targets_verical_configs(targ_idx);

        [fname, fname_base, fname_samples, fname_samples_dat] = RCT.modelFitting.getModelFitFileName(data.subject, data.session, modelSetup);


        fitMAPonly = false;
        if(exist(fname, "File"))
            if(~overwrite)
                mapFile = load(fname,  "summary",  "params_ex", "modelSetup", "fname_samples", "params_map");
                if(false && (~isfield(mapFile, "params_map") || isempty(mapFile.params_map)))
                    fprintf("File found for: %s. Fitting missing MAP estimate...\n", fname_base);
                    fitMAPonly = true;
                    params_ex = mapFile.params_ex;
                    modelSetup = mapFile.modelSetup;
                    fname_samples = mapFile.fname_samples;
                    summary = mapFile.summary;
                else
                    fprintf("File found for: %s. Continuing...\n", fname_base);
                    continue;
                end
            else
                fprintf("File found for: %s. Will overwrite!\n", fname_base);
            end
                
        end

        
        %% get baseline rates for each area on training trials (or test if modelSetup.predictiveBaseline)
        if(~modelSetup.BASELINE_model && modelSetup.baselineCorrectionType ~= "none")
            baselineRates = RCT.modelBuilder.loadBaselineRates(modelSetup);
        else
            baselineRates = []; %
        end

        fprintf("Fitting model: %s.\n", fname_base);

        %%
        if(~fitMAPonly)
            possible_locations = ["FEF" "LIP" "SC"];
            [~,location_idx] = ismember(modelSetup.location, possible_locations);

            coupling_idxs = 0;
            if(modelSetup.includeCouplingInter)
                cs = all_locations(all_locations ~= modelSetup.location);

                if(numel(cs) == 2)
                    coupling_idxs = 3;
                else
                    cs2 = possible_locations(possible_locations ~= modelSetup.location);
                    [~,coupling_idxs] = ismember(cs, cs2);
                end
            end

            rng_samples = str2double(extractBetween(data.session,4,8));
            rng_samples = rng_samples + location_idx * 1e5;
            rng_samples = rng_samples + coupling_idxs * 1e6;
            rng_samples = rng_samples + (modelSetup.includeCouplingLocal + modelSetup.includeTrialwiseLatent*2)* 1e7;
            rng_samples = rng_samples + (modelSetup.runNumber - 1) * 1e8;
    
            modelSetup.rng_seed  = rng_samples;
            modelSetup.GPUs_used = GPUs;
            modelSetup.GPUinfo   = gpuDeviceTable();
        else
            modelSetup.GPUs_used_map = GPUs;
            modelSetup.GPUinfo_map   = gpuDeviceTable();
        end

        %% build gmlm
        if(~fitMAPonly)
            modelSetup.delta_t_sec = data.bin_size_ms * 1e-3;
            modelSetup.llType = [];
        end
        if(~isfield(modelSetup, "llType") || isempty(modelSetup.llType))
            if(data.bin_size_ms > 1)
                modelSetup.llType = "poissExp";
                fprintf("Using Poisson with exponential nonlinearity: bin size = %.1f ms.\n", data.bin_size_ms);
            else
                fprintf("Using truncated Poisson with exponential nonlinearity.\n");
                modelSetup.llType = "truncatedPoissExp";
            end
        end

        [GMLMstructure, trials, modelSetup.spkHistPrior_setup, modelSetup.prior_setup, modelSetup.trialsUsed, modelSetup.couplingShuffle, modelSetup.normalizationSettings] = RCT.modelBuilder.constructGMLMdata14(data, modelSetup, "baselineRates", baselineRates);

        gmlm = kgmlm.GMLM(GMLMstructure, trials, modelSetup.delta_t_sec, modelSetup.llType);
        gmlm.toGPU(GPUs, "useDoublePrecision", modelSetup.doublePrecision);


        %% run HMC
        if(~fitMAPonly)
            HMC_settings = gmlm.setupHMCparams(nWarmup, nSamples, modelSetup.DEBUG);
            HMC_settings.delete_temp_file = true;
            HMC_settings.delete_samples_file = true;
            HMC_settings.sample_H   = false;
            HMC_settings.sample_M   = false;
            HMC_settings.stepSize.maxSteps(1)    = hmc_steps_0;
            HMC_settings.stepSize.maxSteps(2:end)    = hmc_steps;
            HMC_settings.stepSize.delta  = 0.8 * ones(size(HMC_settings.stepSize.schedule, 1),1);

    
            NS_init = [10;10;20;20];
            HMC_settings.stepSize.max_step_size = [logspace(-6,-3,numel(NS_init))'; 1e-2];
            bb = circshift(NS_init,1);
            bb(1) = 1;
            HMC_settings.stepSize.e_init = 1e-6;
            HMC_settings.stepSize.schedule = cat(1, [cumsum(bb) cumsum(NS_init)], HMC_settings.stepSize.schedule);
            HMC_settings.stepSize.schedule(numel(NS_init) + 1,1) = sum(NS_init)+1;

            HMC_settings.stepSize.delta = [0.99*ones(numel(NS_init),1); HMC_delta_0; HMC_delta];
    
            HMC_settings.stepSize.stepL  = HMCstep_L;
    
            tmp_folder = "TempData";
            if(~isfolder(tmp_folder))
                mkdir(tmp_folder);
            end
            HMC_settings.trialLLfile = sprintf("%s/tmp_%s.dat", tmp_folder, fname_base);
            HMC_settings.samplesFile = fname_samples_dat;
    
            rng(rng_samples + 1423);
            params_init = RCT.modelFitting.generateInitialParameters(gmlm, modelSetup);
            
    
    
            HMC_settings.M_init = RCT.modelBuilder.getMInit(params_init, modelSetup, gmlm.dim_M);
    
            if(showPlots)
                figNum = 10;
            else
                figNum = nan;
            end
            opts = gmlm.getComputeOptionsStruct(true);
            res = gmlm.computeLogPosterior(params_init, opts);
            for ii = 1:2
                res = gmlm.computeLogPosterior(params_init, opts, res);
            end
            grad_runs = 100;
            time_0 = tic;
            for ii = 1:grad_runs
                res = gmlm.computeLogPosterior(params_init, opts, res);
            end
            time_final = toc(time_0) * 1e3;
            time_final = time_final / grad_runs;
            fprintf("Timing check: %.2f ms per complete gradient evaluation\n", time_final);
    
            %%
            rng(rng_samples);
            modelSetup.rng_complete_state = rng();
            if(saveDuring)
                HMC_settings.savePartialProgressValid = true;
                HMC_settings.savePartialProgressFile = [sprintf("%s_part1.mat", fname_samples), sprintf("%s_part2.mat", fname_samples)];
                HMC_settings.savePartialProgressN = 100;
    
                printFunc = @RCT.utils.printRegularizedParameterStatus;
                [samples, samples_file_format, summary, HMC_settings] = gmlm.runHMC_simpleLowerRAM( params_init, HMC_settings, "figure", figNum, "printFunc", printFunc, ...
                            "saveSinglePrecision", modelSetup.saveSinglePrecision, "modelInfo", modelSetup); %params_final
                
            else
                [samples, summary, HMC_settings] = gmlm.runHMC_simple( params_init, HMC_settings, "figure", figNum, "saveSinglePrecision", modelSetup.saveSinglePrecision);
                samples_file_format = [];
            end
    
            %% save samples
            save(fname_samples, "-v7.3", "samples", "samples_file_format", "fname_samples_dat",  "params_init",  "HMC_settings", "saveDuring");
       

            %% get parameter fit
    
            % get params of highest post sample (NOT A MAP ESTIMATE!!!!! Just a 'good' sample!)
            [~,ii] = max(samples.log_post((HMC_settings.nWarmup+1):end));
            sample_idx = ii + HMC_settings.nWarmup;
    
            if(saveDuring)
                params_ex = RCT.resultsHandlers.loadSample(fname_samples, sample_idx);
            else
                params_ex = RCT.resultsHandlers.loadSample(samples, sample_idx);
            end
        end

        %% fit a MAP estimate for "good" measure: need to implement a variational approach for the latent variables

        paramStruct_map = [];
        resultsStruct_map = [];
%         fprintf("fitting MAP estimate with hyperparameters set by a posterior save. \n");
%         [paramStruct_map, resultsStruct_map] = gmlm.computeMAP(params_ex, "max_quasinewton_steps", 2000, "max_iters", 2);

        %% save results summary
        save(fname, "-v7.3", "summary",  "params_ex", "modelSetup", "fname_samples", "paramStruct_map", "resultsStruct_map"); %, "paramStruct_map" , "resultsStruct_map", "params_final"

        %% free GPU
        gmlm.freeGPU();

        %% delete partial progress file

        if(isfield(HMC_settings, "savePartialProgressFile")  && isfield(HMC_settings, "savePartialProgressN") && HMC_settings.savePartialProgressN > 0)
            for ff = 1:numel(HMC_settings.savePartialProgressFile)
                if(exist(HMC_settings.savePartialProgressFile(ff), "file"))
                    delete(HMC_settings.savePartialProgressFile(ff));
                end
            end
        end

        %% if baseline model - gets the rates for the training and test trials
        if(modelSetup.BASELINE_model)
            clear gmlm;

            pruneRate = 10;
            sample_idxs = HMC_settings.nWarmup + (1:pruneRate:HMC_settings.nSamples);
            modelSetup_baseline = RCT.modelBuilder.getBaselineModelSetup(modelSetup);
            fname_baseline = RCT.modelFitting.getBaselineModelRateFile(data.subject, data.session, modelSetup_baseline);
            meanRates = RCT.modelFitting.getMeanRatesFromBaseline(data, modelSetup, sample_idxs);
            save(fname_baseline, "meanRates", "sample_idxs", "modelSetup_baseline", "-v7.3");
        end
    end
end
