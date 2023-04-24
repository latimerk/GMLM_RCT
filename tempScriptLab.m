
% clear SC2;
% FFs = ["Results/GMLM/Coupling_none/Mingus/QUARTER123_RCT_GMLM_Mingus_20201211_SC_targV_Rs_S32_R32_C32_iC16_run1.mat"];
% 
% subfieldNames = ["NC"];
% 
% for modelCtr = 1:numel(FFs)
%     load(FFs(modelCtr))
%     RCT.resultsHandlers.setupModelForEvaluation2;
%     RCT.resultsHandlers.allignTotalLL2;
%     SC2.(subfieldNames(modelCtr)) = CS;
% 
%     clear modelSetup
%     clear modelInfo
%     clear HMC_settings
%     clear sample_idx
% end
% 
% save("SC2.mat", "SC2", "-v7.3");


% clear FEF2;
% FFs = ["Results/GMLM/Coupling_none/Mingus/QUARTER123_RCT_GMLM_Mingus_20201211_FEF_targV_Rs_S32_R32_C32_iC16_run1.mat"];
% 
% subfieldNames = ["NC"];
% 
% for modelCtr = 1:numel(FFs)
%     load(FFs(modelCtr))
%     RCT.resultsHandlers.setupModelForEvaluation2;
%     RCT.resultsHandlers.allignTotalLL2;
%     FEF2.(subfieldNames(modelCtr)) = CS;
% 
%     clear modelSetup
%     clear modelInfo
%     clear HMC_settings
%     clear sample_idx
% end
% 
% save("FEF2.mat", "FEF2", "-v7.3");

trainSet = false;

clear LIP2;
FFs = ["Results/GMLM/Coupling_none/Mingus/QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat";
    "Results/GMLM/Coupling_none/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat_part.mat";
    "Results/GMLM/Coupling_none/Mingus/NULL_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat";
    "Results/GMLM/Coupling_none/Mingus/SAMPLES_NULL_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat_part.mat"];

subfieldNames = ["NC", "NC", "NULL", "NULL"];

for modelCtr = 1:numel(FFs)
    if(exist(FFs(modelCtr), "file"))
        try
            load(FFs(modelCtr))
            RCT.resultsHandlers.setupModelForEvaluation2;
            RCT.resultsHandlers.allignTotalLL2;
            LIP2.(subfieldNames(modelCtr)) = CS;
        
            clear modelSetup
            clear modelInfo
            clear HMC_settings
            clear sample_idx
            save("LIP2.mat", "LIP2", "-v7.3");
        catch
            fprintf("Error for file %s\n", FFs(modelCtr));
        end
    end
end

trainSet = true;

% clear TrLIP2;
% FFs = ["Results/GMLM/Coupling_none/Mingus/QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat";
%     "Results/GMLM/Coupling_none/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat_part.mat";
%     "Results/GMLM/Coupling_none/Mingus/NULL_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat";
%     "Results/GMLM/Coupling_none/Mingus/SAMPLES_NULL_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_run1.mat_part.mat"];
% 
% subfieldNames = ["NC", "NC", "NULL", "NULL"];
% 
% for modelCtr = 1:numel(FFs)
%     if(exist(FFs(modelCtr), "file"))
%         try
%             load(FFs(modelCtr))
%             RCT.resultsHandlers.setupModelForEvaluation2;
%             RCT.resultsHandlers.allignTotalLL2;
%             TrLIP2.(subfieldNames(modelCtr)) = CS;
%         
%             clear modelSetup
%             clear modelInfo
%             clear HMC_settings
%             clear sample_idx
%             save("TrLIP2.mat", "TrLIP2", "-v7.3");
%         catch
%             fprintf("Error for file %s\n", FFs(modelCtr));
%         end
%     end
% end