function [baselineRates] = loadBaselineRates(modelSetup_0)

if(modelSetup_0.BASELINE_model)
    baselineRates = [];
else
    if(modelSetup_0.includeCouplingInter)
        locations = [modelSetup_0.location modelSetup_0.interAreaCoupling_locations];
    else
        locations = modelSetup_0.location;
    end
    
    modelSetup_baseline = RCT.modelBuilder.getBaselineModelSetup(modelSetup_0);
            
    for jj = 1:numel(locations)
        modelSetup_baseline.location = locations(jj);
    
        fname_baseline = RCT.modelFitting.getBaselineModelRateFile(modelSetup_0.subject, modelSetup_0.session, modelSetup_baseline);
        if(~exist(fname_baseline, "file"))
            error("No baseline rates found for spike history or coupling!");
        end
        baselineRates.(modelSetup_baseline.location) = load(fname_baseline, "meanRates");
    end
end