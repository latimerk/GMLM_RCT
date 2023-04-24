
addpath ../GMLM_RCT/
addpath ../GMLM
addpath ../GMLM/example/
if(~exist("data", "var"))
%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20200402", "FEF", [], 1);
%     [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20200402", ["FEF" "LIP"]);
    [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", ["FEF"]);
end

dim_R = 16;
shuffle = false;
if(numel(data.trials) > 350)
    %%
%     if(isfield(data.trials(1).Y, "LIP") && shuffle)
%         fprintf("shuffling LIP trials...\n");
%         data_0 = data;
%         order = randperm(numel(data.trials));
%         for ii = 1:numel(data.trials)
%             data.trials(order(ii)).Y.("LIP") = data_0.trials(ii).Y.("LIP");
%         end
%     elseif(isfield(data.trials(1).Y, "LIP") )
%         fprintf("NOT SHUFFLING LIP\n");
%     end

    %%
    data.trials = data.trials(1:8:end);
end


%% get model setup information
rescale = true;
modelSetup.bases = RCT.modelBuilder.setupBasis(data.bin_size_ms, true, rescale);
[modelSetup.stimConfig, modelSetup.responseConfig] = RCT.modelBuilder.getTaskModelSetup(data.TaskInfo);
modelSetup.Ranks= RCT.modelFitting.getDefaultRanks();
modelSetup.includeCouplingLocal   = false;
modelSetup.includeCouplingInter   = false;
modelSetup.includeTrialwiseLatent = true;
modelSetup.includeTrialwiseLatent_allStimEvents = true;
modelSetup.includeTrialwiseLatent_scaleFunc = @(dim_M) 1 / sqrt(dim_M);
modelSetup.constant_scale = 0.1;

modelSetup.reprojectLatents    = false;%runNum == 2;
modelSetup.removeUnusedLatents = false;

modelSetup.includeCouplingInter_trialDim = false; 
modelSetup.includeCouplingLocal_trialDim = false;

modelSetup.Ranks.stimulus = dim_R;
modelSetup.Ranks.response = dim_R;
modelSetup.Ranks.coupling_local = dim_R;
modelSetup.Ranks.coupling_inter = dim_R;

modelSetup.location = "FEF";
modelSetup.targets_verical = true;
modelSetup.fitHalf = 0;
modelSetup.useGibbsStep = false;
modelSetup.delta_t_sec = data.bin_size_ms * 1e-3;

[GMLMstructure, trials, prior_setup, trialsUsed, latentBasis] = RCT.modelBuilder.constructGMLMdata9(data, modelSetup);

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
if(~exist("GPUs", "var"))
    if(gpuDeviceCount() > 2)
        
        if(exist("runNum", "var"))
            if(runNum == 1)
                GPUs = [0 1];
            elseif(runNum == 2)
                GPUs = [2 3];
            else
                GPUs = [0 1 2];
            end
        else
            GPUs = [2 3];
        end
    elseif(gpuDeviceCount() > 0)
        GPUs = 0;
    else
        error("no GPUs");
    end
end


if(~gmlm.isOnGPU())
    gmlm.toGPU(GPUs, "useDoublePrecision", true);
end

for jj = 1:gmlm.dim_J
    gmlm.setDimR(jj,dim_R);
end
%%
params = gmlm.getRandomParamStruct();
opts = gmlm.getComputeOptionsStruct(true);
res = gmlm.getEmptyResultsStruct(opts);
tic;
for ii = 1:100
    res = gmlm.computeLogPosterior(params, opts, res);
end
toc


    nWarmup = 40e3;
    nSamples = 80e3;
HMC_settings = gmlm.setupHMCparams(nWarmup, nSamples);
HMC_settings.sample_M   = false;
HMC_settings.sample_H   = false;
HMC_settings.fitMAP= [];
HMC_settings.stepSize.maxSteps(:) = 500;
HMC_settings.delete_temp_file = true;
HMC_settings.delete_samples_file = true;
HMC_settings.stepSize.scaleDuringWarmup = false;
HMC_settings.sample_M_setScale = false;
HMC_settings.samplesFile = "TempData/SAMPLES.mat";
HMC_settings.M_init = RCT.modelBuilder.getMInit(params, modelSetup, gmlm.dim_M);

if(runNum == 1)
    HMC_settings.M_est.first_sample  = []; 
    HMC_settings.M_est.samples       = [];
    
    HMC_settings.stepSize.schedule   = [2      20000;
                                        20001  35000]; %each row gives a range of trials to estimate step size (restarts estimation at each sample = schedule(ii,1))

    HMC_settings.stepSize.scaleRanges = [35001 nSamples];
end

if(runNum == 1)
    %%
    modelType = "plain latent: fixed M";
   
    HMC_settings.trialLLfile = "TempData/s1.mat";
    [samples1, summary1] = gmlm.runHMC_simple( params, HMC_settings, "figure", 11);

    samplesPart1.log_like = samples1.log_like;
    samplesPart1.log_post = samples1.log_post;
    samplesPart1.W = samples1.W;
    save("samples1.mat", "samples1", "summary1", "modelType", "samplesPart1", "-v7.3");
elseif(runNum == 2)
    modelType = "plain latent: fit M";
    HMC_settings.trialLLfile = "TempData/s2.mat";
    [samples2, summary2] = gmlm.runHMC_simple( params, HMC_settings, "figure", 12);

    samplesPart2.log_like = samples2.log_like;
    samplesPart2.log_post = samples2.log_post;
    samplesPart2.W        = samples2.W;
    save("samples2.mat", "samples2", "summary2", "modelType", "samplesPart2", "-v7.3");
elseif(runNum == 3)
    %%
    modelType = "plain latent 2";
   
    HMC_settings.trialLLfile = "TempData/s3.mat";
    [samples3, summary3] = gmlm.runHMC_simple( params, HMC_settings, "figure", 11);

    samplesPart3.log_like = samples3.log_like;
    samplesPart3.log_post = samples3.log_post;
    samplesPart3.W = samples3.W;
    save("samples3.mat", "samples3", "summary3", "modelType", "samplesPart3", "-v7.3");
end



return;
%%
HMC_settings.samplesFile = "TempData/s2.mat";
HMC_settings.stepSize.maxSteps = 200;
HMC_settings.sample_M   = false;
HMC_settings.sample_H   = false;
%%

% params = gmlm.getRandomParamStruct();
% opts = gmlm.getComputeOptionsStruct(true, "trial_weights", true);
% 
% bb = nan(102,1);
% for kk = 1:102
%     K = kk;
%     fprintf("K = %d\n", K);
%     opts.trial_weights(:) = 0;
%     opts.trial_weights(randperm(gmlm.dim_M, K)) = gmlm.dim_M/K;
%     % opts.Groups(3).dV = false;
%     % opts.Groups(3).dT(1) = false;
%     % opts.Groups(3).dT(2) = false;
%     
%     res = gmlm.getEmptyResultsStruct(opts);
%     aa = tic;
%     for ii = 1:100
%         res = gmlm.computeLogPosterior(params, opts, res);
%     end
%     bb(kk) = toc(aa);
% end
