
% RCT.utils.setupGMLMpaths()
% quickRun
% [~,summary,~,paramStruct] = gmlm.runHMC_simple( params, HMC_settings, "figure", 2, "printFunc", printFunc, "runWAIC", false);

% save("quickRunInit.mat", "summary", "paramStruct", "-v7.3");
quickRun
load("quickRunInit.mat", "summary", "paramStruct");
%%
e_0 = min(0.05, summary.HMC_state.stepSize.e_bar);

HMC_settings = gmlm.setupHMCparams(0, 11e3);
HMC_settings.stepSize.schedule = [];
HMC_settings.stepSize.max_step_size = 1;
HMC_settings.stepSize.stepL = 100;

HMC_settings.delete_temp_file = true;
HMC_settings.delete_samples_file = true;
HMC_settings.samplesFile = "TempData/SAMPLES_quickRun.mat";
HMC_settings.trialLLfile = "TempData/quickRun.mat";
HMC_settings.M_init = RCT.modelBuilder.getMInit(params, modelSetup, gmlm.dim_M);


Nsteps = 60:20:200;

for ss = 1:numel(Nsteps)
    HMC_settings.stepSize.maxSteps = Nsteps(ss);
    HMC_settings.stepSize.e_init = e_0;
    HMC_settings.stepSize.e_0 = HMC_settings.stepSize.e_init;
   
    samples_test(ss) = gmlm.runHMC_simple( paramStruct, HMC_settings, "figure", 10 + ss, "printFunc", printFunc, "runWAIC", false);
end
save("quickRunResults.mat", "-v7.3", "samples_test");