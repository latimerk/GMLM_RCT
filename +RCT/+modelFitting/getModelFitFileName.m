function [fname, fname_base, fname_samples, fname_samples_dat, baseFolder] = getModelFitFileName(subject, session_name, modelSetup)

[folders, subject, ~, session_name, valid_modelSetup.locations] = RCT.dataHandlers.checkSession([], subject, session_name);
if(~isfolder(folders.data.processed.(subject)))
    mkdir(folders.data.processed.(subject));
end

if(ischar(modelSetup.location))
    modelSetup.location = string(modelSetup.location);
end
if(numel(modelSetup.location) > 1 ||  ~isstring(modelSetup.location))
    error("Invalid recording modelSetup.location specified");
end
if(~ismember(modelSetup.location, valid_modelSetup.locations))
    error("Invalid recording modelSetup.location specified");
end

if(~isempty(modelSetup.interAreaCoupling_locations))
    if(ischar(modelSetup.interAreaCoupling_locations))
        modelSetup.interAreaCoupling_locations = string(modelSetup.interAreaCoupling_locations);
    elseif(~isstring(modelSetup.interAreaCoupling_locations))
        error("Invalid coupling modelSetup.location specified");
    end

    if(ismember(modelSetup.location, modelSetup.interAreaCoupling_locations))
        error("interarea modelSetup.locations cannot include current recording area");
    elseif(~all(ismember(modelSetup.interAreaCoupling_locations, valid_modelSetup.locations)))
        error("Invalid coupling modelSetup.location specified");
    end
end


if(isfield(modelSetup, "BASELINE_model") && modelSetup.BASELINE_model)
    baseFolder = folders.results.GMLM.baseline.(subject);
elseif(~modelSetup.includeCouplingLocal && isempty(modelSetup.interAreaCoupling_locations))
    baseFolder = folders.results.GMLM.noCoupling.(subject);
elseif(modelSetup.includeCouplingLocal && isempty(modelSetup.interAreaCoupling_locations))
    baseFolder = folders.results.GMLM.localCoupling.(subject);
elseif( ~isempty(modelSetup.interAreaCoupling_locations))
    baseFolder = folders.results.GMLM.allCoupling.(subject);
end
baseFolder = sprintf("%s/S_%s/", baseFolder, session_name);

if(~isfolder(baseFolder))
    mkdir(baseFolder);
end

if(modelSetup.targets_verical)
    targets_str = "_targV";
else
    targets_str = "_targH";
end

localCoupling_str = "";
if(modelSetup.includeCouplingLocal)
    localCoupling_str = "_LC";
end

trialwise_str = "";
if(modelSetup.includeTrialwiseLatent)
    if(modelSetup.includeTrialwiseLatent_allStimEvents)
        trialwise_str = "_TrLsA";
    else
        trialwise_str = "_TrLs1";
    end
end

interAreaCoupling_str = "";
if(~isempty(modelSetup.interAreaCoupling_locations))
    if(isfield(modelSetup, "interAreaCoupling_shuffle") && modelSetup.interAreaCoupling_shuffle)
        if(modelSetup.interAreaCoupling_shuffle == 2)
            interAreaCoupling_str = "_ShuffledSacciC";
        else
            interAreaCoupling_str = "_ShuffledStimiC";
        end
    else
        interAreaCoupling_str = "_iC";
    end
    for ii = 1:numel(modelSetup.interAreaCoupling_locations)
        interAreaCoupling_str = sprintf("%s_%s", interAreaCoupling_str, modelSetup.interAreaCoupling_locations(ii));
    end
end

DEBUG_str = "";
if(modelSetup.DEBUG)
    DEBUG_str = "DEBUG_";
end

REGULARIZER_str = "";
if(isfield(modelSetup, "useRidge") && modelSetup.useRidge)
    REGULARIZER_str = "RIDGE_";
end


half_str= "";
if(isfield(modelSetup, "fitQuarter") && all(modelSetup.fitQuarter > 0))
    qs = unique(modelSetup.fitQuarter);
    half_str= "QUARTER";
    for ii = 1:numel(qs)
        half_str= sprintf("%s%d", half_str, qs(ii));
    end
    half_str= sprintf("%s_", half_str);
elseif(isfield(modelSetup, "fitHalf") && modelSetup.fitHalf > 0)
    half_str= sprintf("HALF%d_", modelSetup.fitHalf);
end

prec_str = "";
if(~modelSetup.doublePrecision)
    prec_str = "SINGLEPRECISION_";
end

NULL_str = "";
if(isfield(modelSetup, "NULL_Model") && modelSetup.NULL_Model)
    NULL_str = "NULL_";
elseif(isfield(modelSetup, "BASELINE_model") && modelSetup.BASELINE_model)
    NULL_str = "BASELINE_";
elseif(isfield(modelSetup, "baselineCorrectionType") && lower(modelSetup.baselineCorrectionType) ~= "none")
    NULL_str = sprintf("BC%s_", lower(modelSetup.baselineCorrectionType));
end

run_str = sprintf("_run%d", modelSetup.runNumber);

Rank_str = sprintf("_Rs_S%d_R%d_C%d_iC%d", modelSetup.Ranks.stimulus, modelSetup.Ranks.response, modelSetup.Ranks.coupling_local, modelSetup.Ranks.coupling_inter);

base_str = sprintf("%s%s%s%s%sRCT_GMLM_%s_%s_%s", NULL_str, REGULARIZER_str, DEBUG_str, half_str, prec_str, subject, session_name, modelSetup.location);

fname_base = sprintf("%s%s%s%s%s%s%s", base_str, targets_str, Rank_str, trialwise_str, localCoupling_str, interAreaCoupling_str, run_str);
fname = sprintf("%s/%s.mat", baseFolder, fname_base); 
fname_samples = sprintf("%s/SAMPLES_%s.mat", baseFolder, fname_base); 
fname_samples_dat = sprintf("%s/SAMPLESDAT_%s.dat", baseFolder, fname_base); 