
addpath ../GMLM_RCT/
addpath ../GMLM_dmc/
addpath ../GMLM_dmc/example/
clear modelSetup;

modelSetup.location = "FEF";
modelSetup.targets_verical = true;

removeInvalidCells = true; % will need to setup a thing for each dataset to remove any obviously missorted cells

downsample = 5; %ms bins 
includeSpkHistory = false;

GPUs = 0; % which GPUs to use

[data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", modelSetup.location,removeInvalidCells, downsample );

%% get model setup information
[modelSetup.stimConfig, modelSetup.responseConfig] = RCT.modelBuilder.getTaskModelSetup(data.TaskInfo);

modelSetup.bases = RCT.modelBuilder.setupBasis();
modelSetup.Ranks = RCT.modelFitting.getDefaultRanks();
modelSetup.kernelTimescales =  RCT.modelBuilder.getDefaultKernelTimescales();

if(~includeSpkHistory)
    modelSetup.bases.spkHist.B = [];
end

modelSetup.BASELINE_model = false;
modelSetup.flatPrior = true;

modelSetup.constant_scale = 0.1;

modelSetup.includeCouplingLocal   = false;
modelSetup.includeCouplingInter   = false;
modelSetup.includeTrialwiseLatent = false;
modelSetup.includeTrialwiseLatent_allStimEvents = false;
modelSetup.includeTrialwiseLatent_log_scaleFunc = @(dim_M) -1/2*log(dim_M);
modelSetup.interAreaCoupling_shuffle = 0;

modelSetup.baselineCorrectionType = "none";
modelSetup.baselineCorrectSpikeHistory = false;


modelSetup.reprojectLatents    = false;
modelSetup.removeUnusedLatents = false;

modelSetup.includeCouplingInter_trialDim = false; 
modelSetup.includeCouplingLocal_trialDim = false;

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

%%

if(gpuDeviceCount() < 1)
    error("no GPUs");
end
if(~gmlm.isOnGPU())
    gmlm.toGPU(GPUs, "useDoublePrecision", modelSetup.doublePrecision);
end

%%
numWarmup = 10e3;
numSamples = 20e3;

params = gmlm.getRandomParamStruct();
HMC_settings = gmlm.setupHMCparams(numWarmup, numSamples);
HMC_settings.samplesFile = "TempData/SAMPLES_quickRun.mat";
HMC_settings.trialLLfile = "TempData/quickRun.mat";
HMC_settings.M_init = RCT.modelBuilder.getMInit(params, modelSetup, gmlm.dim_M);
printFunc = @RCT.utils.printRegularizedParameterStatus;

[samples1, summary1] = gmlm.runHMC_simple( params, HMC_settings, "figure", 11, "printFunc", printFunc, "runWAIC", false);