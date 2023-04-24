function [fname, baseFolder] = getBaselineModelRateFile(subject, session_name, modelSetup_baseline)

 [~, fname_base, ~, ~, baseFolder]= RCT.modelFitting.getModelFitFileName(subject, session_name, modelSetup_baseline);

 baseFolder = sprintf("%s/baselineRates/", baseFolder);
 if(~isfolder(baseFolder))
     mkdir(baseFolder);
 end

 fname = sprintf("%s/%s_rates.mat", baseFolder, fname_base);