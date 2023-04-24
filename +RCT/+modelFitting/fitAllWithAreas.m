function [] = fitAllWithAreas(locations,  GPUs, varargin)


p = inputParser;
p.CaseSensitive = false;

addRequired(p, "locations",     @(A) (ischar(A) || isstring(A)) && ~isempty(A));
addRequired(p, "GPUs",     @(A) isnumeric(A)  && ~isempty(A));
p = RCT.modelFitting.setupDefaultModelFittingArgs(p);

parse(p, locations, GPUs, varargin{:});
% then set/get all the inputs out of this structure
GPUs  = p.Results.GPUs;
locations  = p.Results.locations;
Ranks                    = p.Results.Ranks;
includeCouplingLocal     = p.Results.includeCouplingLocal;
includeCouplingInter     = numel(locations) > 1;
includeTrialwiseLatent     = p.Results.includeTrialwiseLatent;
includeTrialwiseLatent_allStimEvents     = p.Results.includeTrialwiseLatent_allStimEvents;
includeTrialwiseLatent_log_scaleFunc = p.Results.includeTrialwiseLatent_log_scaleFunc;

fitHalf                = p.Results.fitHalf;
fitQuarter                = p.Results.fitQuarter;

showPlots                = p.Results.showPlots;
saveDuring = p.Results.saveDuring;
runNumber                = p.Results.runNumber;
[~,IT]  = unique(p.Results.targets_verical_configs);
targets_verical_configs = p.Results.targets_verical_configs(sort(IT));
overwrite                = p.Results.overwrite;
subjects                = p.Results.subjects;
DEBUG                = p.Results.DEBUG;
saveSinglePrecision                = p.Results.saveSinglePrecision;
doublePrecision                = p.Results.doublePrecision;
hmc_steps = p.Results.maxHMCsteps;
hmc_steps_0 = p.Results.maxHMCsteps_0;
HMCstep_L = p.Results.HMCstep_L;
HMC_delta = p.Results.HMC_delta;
HMC_delta_0 = p.Results.HMC_delta_0;
nSamples = p.Results.nSamples;
nWarmup = p.Results.nWarmup;
useGibbsStep = p.Results.useGibbsStep;
downsample = p.Results.downsample;
shuffleCoupling = p.Results.shuffleCoupling;
useRidge = p.Results.useRidge;
useNull = p.Results.useNull;
useBaseline = p.Results.useBaseline;


baselineCorrectionType = p.Results.baselineCorrectionType;

folders = RCT.dataHandlers.getFolders();
if(ischar(locations))
    locations = string(locations);
end
if(ischar(subjects))
    subjects = string(subjects);
end
if(isempty(subjects))
    subjects = folders.subjects;
end

if(~all(ismember(subjects, folders.subjects)))
    error("Invalid subjects requested");
end

%% for each session containing all the requested areas
for sub_idx = 1:numel(subjects)
    subject_name = subjects(sub_idx);

    % for each session
    for sess_idx = 1:numel(folders.sessions.(subject_name))
        if(all(ismember(locations, folders.sessions.(subject_name)(sess_idx).locations)))
            data = RCT.dataHandlers.loadData(subject_name, folders.sessions.(subject_name)(sess_idx).name, locations, [], downsample);

            RCT.modelFitting.runModelFitting(data, GPUs,  locations(1), "Ranks", Ranks, ...
                "includeTrialwiseLatent", includeTrialwiseLatent, "includeTrialwiseLatent_allStimEvents", includeTrialwiseLatent_allStimEvents, "includeTrialwiseLatent_log_scaleFunc", includeTrialwiseLatent_log_scaleFunc,  ...
                "includeCouplingLocal", includeCouplingLocal, "includeCouplingInter", includeCouplingInter, ...
                "useGibbsStep", useGibbsStep, "useRidge", useRidge, "useNull", useNull, "useBaseline", useBaseline, "baselineCorrectionType", baselineCorrectionType,  ...
                "targets_verical_configs", targets_verical_configs, ...
                "showPlots", showPlots, "saveDuring", saveDuring, "runNumber", runNumber, "overwrite", overwrite, "debug", DEBUG, "doublePrecision", doublePrecision, "saveSinglePrecision", saveSinglePrecision, "fitHalf", fitHalf, "fitQuarter", fitQuarter, ...
                "maxHMCsteps", hmc_steps, "maxHMCsteps_0", hmc_steps_0, "HMCstep_L", HMCstep_L, "HMC_delta", HMC_delta, "HMC_delta_0", HMC_delta_0, "nSamples", nSamples, "nWarmup", nWarmup, "shuffleCoupling", shuffleCoupling);
        end
    end
end
