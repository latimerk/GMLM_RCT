% clear SC;
% FFs = ["Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_SC_targV_Rs_S32_R32_C32_iC16_LC_iC_FEF_LIP_run1.mat_part.mat";
%        "Results/GMLM/Coupling_local/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_SC_targV_Rs_S32_R32_C32_iC16_LC_run1.mat_part.mat";
%        "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_SC_targV_Rs_S32_R32_C32_iC16_LC_iC_FEF_run1.mat_part.mat";
%        "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_SC_targV_Rs_S32_R32_C32_iC16_LC_iC_LIP_run1.mat_part.mat"];
% 
% subfieldNames = ["FC", "LC", "PC_FEF", "PC_LIP"];
% 
% for modelCtr  = 1:numel(FFs)
%     if(exist(FFs(modelCtr), "file"))
%         load(FFs(modelCtr))
%         RCT.resultsHandlers.setupModelForEvaluation2;
%         RCT.resultsHandlers.allignTotalLL2;
%         SC.(subfieldNames(modelCtr)) = CS;
%     
%         clear modelSetup
%         clear modelInfo
%         clear HMC_settings
%         clear sample_idx
%         save("SC.mat", "SC", "-v7.3");
%     end
% 
% end


trainSet = false;

clear LIP;
FFs = ["Results/GMLM/Coupling_local/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_run1.mat_part.mat";
       "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_iC_FEF_SC_run1.mat_part.mat";
       "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_iC_FEF_run1.mat_part.mat";
       "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_iC_SC_run1.mat_part.mat"];

subfieldNames = ["LC", "FC",  "PC_FEF", "PC_LIP"];

for modelCtr  = 1:numel(FFs)
    if(exist(FFs(modelCtr), "file"))
        try
            load(FFs(modelCtr))
            RCT.resultsHandlers.setupModelForEvaluation2;
            RCT.resultsHandlers.allignTotalLL2;
            LIP.(subfieldNames(modelCtr)) = CS;
        
            clear modelSetup
            clear modelInfo
            clear HMC_settings
            clear sample_idx
            save("LIP.mat", "LIP", "-v7.3");
        catch
            fprintf("Error for file %s\n", FFs(modelCtr));
        end
    end
end


% trainSet = true;
% 
% clear TrLIP;
% FFs = ["Results/GMLM/Coupling_local/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_run1.mat_part.mat";
%        "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_iC_FEF_SC_run1.mat_part.mat";
%        "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_iC_FEF_run1.mat_part.mat";
%        "Results/GMLM/Coupling_all/Mingus/SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_iC_SC_run1.mat_part.mat"];
% 
% subfieldNames = ["LC", "FC",  "PC_FEF", "PC_LIP"];
% 
% for modelCtr  = 1:numel(FFs)
%     if(exist(FFs(modelCtr), "file"))
%         try
%             load(FFs(modelCtr))
%             RCT.resultsHandlers.setupModelForEvaluation2;
%             RCT.resultsHandlers.allignTotalLL2;
%             TrLIP.(subfieldNames(modelCtr)) = CS;
%         
%             clear modelSetup
%             clear modelInfo
%             clear HMC_settings
%             clear sample_idx
%             save("TrLIP.mat", "TrLIP", "-v7.3");
%         catch
%             fprintf("Error for file %s\n", FFs(modelCtr));
%         end
%     end
% end