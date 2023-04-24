function [fname] = getRawDataFileName(folders, subject, session_name, location)
[folders, subject, ~, session_name, location] = RCT.dataHandlers.checkSession(folders, subject, session_name, location);

if(~isfolder(folders.data.raw.(subject)))
    mkdir(folders.data.raw.(subject));
end
fname = [sprintf("%s/Analyzed_merged_%s_%s_%s.mat", folders.data.raw.(subject), subject, session_name, location);
         sprintf("%s/merged_%s_%s_%s.mat", folders.data.raw.(subject), subject, session_name, location)];