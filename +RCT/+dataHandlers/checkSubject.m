function [folders, subject] = checkSubject(folders, subject)

%% if folders is not given, gets it
if(isempty(folders))
    folders = RCT.dataHandlers.getFolders();
end

%% checks subject
% converts subject to string
if(isempty(subject) || (isnumeric(subject) && (subject < 1 || subject > numel(folders.subjects))))
    error("Invalid subject.");
elseif(isnumeric(subject))
    subject = folders.subjects(subject);
end
if(ischar(subject))
    subject = string(subject);
end

if(all(folders.subjects ~= subject))
    error("Subject not found.");
end