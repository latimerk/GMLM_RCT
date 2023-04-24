function [fname, fname_curation] = getProcessedDataFileName(folders, subject, session_name)
[folders, subject, ~, session_name] = RCT.dataHandlers.checkSession(folders, subject, session_name);
if(~isfolder(folders.data.processed.(subject)))
    mkdir(folders.data.processed.(subject));
end

fname = sprintf("%s/Session_%s_complete.mat", folders.data.processed.(subject), session_name);
fname_curation = sprintf("%s/Session_%s_complete_curation.mat", folders.data.processed.(subject), session_name);