function [modelSetup_baseline] = getBaselineModelSetup(modelSetup_0)

modelSetup_baseline.targets_verical = modelSetup_0.targets_verical;
modelSetup_baseline.fitQuarter = modelSetup_0.fitQuarter;
modelSetup_baseline.fitHalf = modelSetup_0.fitHalf;
modelSetup_baseline.doublePrecision = modelSetup_0.doublePrecision;
modelSetup_baseline.runNumber = modelSetup_0.runNumber;
modelSetup_baseline.Ranks  = modelSetup_0.Ranks;
modelSetup_baseline.location  = modelSetup_0.location;
modelSetup_baseline.DEBUG  = modelSetup_0.DEBUG;

modelSetup_baseline.BASELINE_model = true;
modelSetup_baseline.interAreaCoupling_locations = [];
modelSetup_baseline.includeCouplingInter = false;
modelSetup_baseline.includeCouplingLocal = false;
modelSetup_baseline.includeTrialwiseLatent = false;
modelSetup_baseline.NULL_model = false;

modelSetup_baseline.stimConfig = modelSetup_0.stimConfig;
modelSetup_baseline.responseConfig = modelSetup_0.responseConfig;

TaskInfo_dummy.categories = [1;1;1;2;2;2];
TaskInfo_dummy.directions = [75;135;195;255;315;15];
[~, ~, modelSetup_baseline.fixationConfig] = RCT.modelBuilder.getTaskBaselineModelSetup(TaskInfo_dummy);