
function [stimConfig, responseConfig, fixationConfig] = getTaskBaselineModelSetup(TaskInfo)

[stimConfig, responseConfig, fixationConfig] = RCT.modelBuilder.getTaskModelSetup(TaskInfo);

end