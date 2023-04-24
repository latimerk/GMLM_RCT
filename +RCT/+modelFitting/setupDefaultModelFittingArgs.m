function [p] = setupDefaultModelFittingArgs(p)


addParameter(p, "Ranks"    ,  RCT.modelFitting.getDefaultRanks(),    @(A) isstruct(A) && isfield(A, "stimulus") && isfield(A, "coupling_local") && isfield(A, "coupling_inter") && isfield(A, "response"));
addParameter(p, "fitHalf",  0,    @(aa) ismember(aa, [0 1 2]));
addParameter(p, "fitQuarter",  [1 2 3],    @(aa) all(ismember(aa, [0 1 2 3 4])));
addParameter(p, "includeCouplingInter",  true,    @islogical);
addParameter(p, "includeCouplingLocal",  true,    @islogical);
addParameter(p, "includeTrialwiseLatent",  true,    @islogical);
addParameter(p, "includeTrialwiseLatent_allStimEvents",  true,    @islogical);
addParameter(p, "includeTrialwiseLatent_log_scaleFunc",  @(dim_M) -0.5*log(dim_M) ,    @(a)isa(a, "function_handle")    ); 
addParameter(p, "showPlots",  false,    @islogical);
addParameter(p, "saveDuring",  true,    @islogical);
addParameter(p, "runNumber",  1,    @(A) isnumeric(A) && isscalar(A));
addParameter(p, "targets_verical_configs",  true,    @islogical);
addParameter(p, "overwrite",  false,    @islogical);
addParameter(p, "subjects", [],    @(A) (ischar(A) || isstring(A))  || ~isempty(A));
addParameter(p, "DEBUG",  false,    @islogical);
addParameter(p, "doublePrecision",  true,    @islogical);
addParameter(p, "saveSinglePrecision",  true,    @islogical);
addParameter(p, "maxHMCsteps",  100,    @(aa) isnumeric(aa) && fix(aa) == aa && aa > 0);
addParameter(p, "maxHMCsteps_0",  100,    @(aa) isnumeric(aa) && fix(aa) == aa && aa > 0);
addParameter(p, "HMCstep_L",  1,    @(aa) isnumeric(aa) && aa > 0);
addParameter(p, "HMC_delta",  0.9,    @(aa) isnumeric(aa) && aa > 0 && aa < 1);
addParameter(p, "HMC_delta_0",  0.8,    @(aa) isnumeric(aa) && aa > 0 && aa < 1);
addParameter(p, "nSamples",  40e3,    @(aa) isnumeric(aa) && aa > 0);
addParameter(p, "nWarmup",  10e3,    @(aa) isnumeric(aa) && aa > 0);
addParameter(p, "downsample",  1,    @isnumeric);
addParameter(p, "useGibbsStep",  false,    @islogical);
addParameter(p, "shuffleCoupling",  false,    @(aa)ismember(aa,[false 1 2]));
addParameter(p, "useRidge",  false,    @islogical);
addParameter(p, "useNull",  false,    @islogical);
addParameter(p, "useBaseline",  false,    @islogical);
addParameter(p, "baselineCorrectionType",  "none",    @(a)ismember(lower(a), ["none", "subtract", "rate", "zscore", "bernoulli"]));