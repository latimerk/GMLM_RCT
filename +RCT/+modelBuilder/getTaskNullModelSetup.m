
function [stimConfig, responseConfig] = getTaskNullModelSetup(TaskInfo)

[stimConfig, responseConfig] = RCT.modelBuilder.getTaskModelSetup(TaskInfo, "includeDirection", false, "includeCategory", false, "includeTargetLocationsOnset", false, ...
            "includeTargetLocationsFixOffset",  false, "includeResponseLocations", false, "includeResponseTargetColor",  false, "includeResponseLocationTargetInteractions",  false);

end