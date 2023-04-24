
addpath ../GMLM_RCT/
addpath ../GMLM
addpath ../GMLM/example/
if(~exist("data", "var"))
%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20200402", "FEF", [], 1);
%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", ["FEF" "LIP" "SC"]);

%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", ["LIP" "FEF" "SC"]);
%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", ["LIP" ]);
%     data.trials = data.trials(1:400);
%     data.trials = data.trials(~ismember(1:400, 4:4:400));

%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", ["LIP" "SC"], true, 2);
%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", ["FEF"], true, 2);
    load("quickRunData.mat", "data", "session_name");

    if(numel(data.trials) > 350)
    
        %%
        data.trials = data.trials(1:12:400);
    end

end

dim_R = 2;

%% get model setup information
clear modelSetup;
modelSetup.bases = RCT.modelBuilder.setupBasis();
modelSetup.BASELINE_model = false;

if(isfield(modelSetup, "BASELINE_model") && modelSetup.BASELINE_model)
    [modelSetup.stimConfig, modelSetup.responseConfig, modelSetup.fixationConfig] = RCT.modelBuilder.getTaskBaselineModelSetup(data.TaskInfo);
    modelSetup.bases.spkHist.B = [];
elseif(isfield(modelSetup, "NULL_model") && modelSetup.NULL_model)
    [modelSetup.stimConfig, modelSetup.responseConfig] = RCT.modelBuilder.getTaskNullModelSetup();
else
    [modelSetup.stimConfig, modelSetup.responseConfig] = RCT.modelBuilder.getTaskModelSetup(data.TaskInfo);
end

modelSetup.flatPrior = true;
modelSetup.Ranks = RCT.modelFitting.getDefaultRanks();%3, 2, 2, 2, 2);

modelSetup.includeCouplingLocal   = false;
modelSetup.includeCouplingInter   = false;
modelSetup.includeTrialwiseLatent = false;
modelSetup.includeTrialwiseLatent_allStimEvents = true;
modelSetup.includeTrialwiseLatent_log_scaleFunc = @(dim_M) -1/2*log(dim_M);
modelSetup.constant_scale = 0.1;
modelSetup.interAreaCoupling_shuffle = 0;

modelSetup.baselineCorrectionType = "none";
modelSetup.baselineCorrectSpikeHistory = false;


modelSetup.kernelTimescales =  RCT.modelBuilder.getDefaultKernelTimescales();

modelSetup.reprojectLatents    = false;
modelSetup.removeUnusedLatents = false;

modelSetup.includeCouplingInter_trialDim = false; 
modelSetup.includeCouplingLocal_trialDim = false;


modelSetup.location = "FEF";
modelSetup.targets_verical = true;
modelSetup.fitHalf = 0;
modelSetup.useGibbsStep = false;
modelSetup.delta_t_sec = data.bin_size_ms * 1e-3;
modelSetup.doublePrecision = true;

modelSetup.fitQuarter = 1:3;
modelSetup.runNumber = 1;
modelSetup.DEBUG = false;
modelSetup.subject = data.subject;
modelSetup.session = data.session;


modelSetup.fitQuarter = 0;
[GMLMstructure, trials, modelSetup.spkHistPrior_setup, modelSetup.prior_setup, trialsUsed, couplingSetup, modelSetup.normalizationSettings] = RCT.modelBuilder.constructGMLMdata14(data, modelSetup);


if(modelSetup.delta_t_sec <= 1e-3)
    fprintf("Using truncated Poisson exp\n");
    ll_type = "truncatedPoissExp";
else
    fprintf("Using Poisson exp (bin size = %.1f ms)\n", data.bin_size_ms);
    ll_type = "poissExp";
end
% ll_type = "poissExp";
% ll_type = "poissSoftRec";
gmlm  = kgmlm.GMLM(GMLMstructure, trials, modelSetup.delta_t_sec, ll_type);

opts = gmlm.getComputeOptionsStruct(true);
params = RCT.modelFitting.generateInitialParameters(gmlm, modelSetup);
res = gmlm.getEmptyResultsStruct(opts);


%%
if(gpuDeviceCount() > 1)
    GPUs = [0 1 2 3];
elseif(gpuDeviceCount() > 0)
    GPUs = 0;
else
    error("no GPUs");
end


if(~gmlm.isOnGPU())
    gmlm.toGPU(GPUs, "useDoublePrecision", modelSetup.doublePrecision);
end

% for jj = 1:gmlm.dim_J
%     gmlm.setDimR(jj, dim_R);
% end

%%
params = RCT.modelFitting.generateInitialParameters(gmlm, modelSetup);
% params = gmlm.getRandomParamStruct();
params.B(:) = randn(gmlm.dim_B, gmlm.dim_P)*0.1;

opts = gmlm.getComputeOptionsStruct(true);
opts_none = gmlm.getComputeOptionsStruct(false);

res = gmlm.getEmptyResultsStruct(opts);
for ii = 1:2
    res = gmlm.computeLogPosterior(params, opts, res);
end
tic;
for ii = 1:100
    res = gmlm.computeLogPosterior(params, opts, res);
end
toc


HMC_settings = gmlm.setupHMCparams(10e3, 5e3);
HMC_settings.stepSize.e_init = 1e-4;
HMC_settings.stepSize.schedule(1) = 10;
HMC_settings.stepSize.scaleRanges = cat(1, [1 HMC_settings.stepSize.schedule(1)], HMC_settings.stepSize.scaleRanges);

HMC_settings.stepSize.maxSteps(:) = 100;
HMC_settings.delete_temp_file = true;
HMC_settings.delete_samples_file = true;
HMC_settings.samplesFile = "TempData/SAMPLES_quickRun.mat";
HMC_settings.trialLLfile = "TempData/quickRun.mat";
HMC_settings.M_init = RCT.modelBuilder.getMInit(params, modelSetup, gmlm.dim_M);
printFunc = @RCT.utils.printRegularizedParameterStatus;

return;
%%
[samples1, summary1] = gmlm.runHMC_simple( params, HMC_settings, "figure", 11, "printFunc", printFunc, "runWAIC", false);
%%

opts = gmlm.getComputeOptionsStruct(true);
opts_none = gmlm.getComputeOptionsStruct(false);



tic;
for ii = 1:1000
    res = gmlm.computeLogPosterior(params, opts, res);
end
toc
tic;
for ii = 1:1000
    res = gmlm.computeLogPosterior(params, opts_none, res);
end
toc