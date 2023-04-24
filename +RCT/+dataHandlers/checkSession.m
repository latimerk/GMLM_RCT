function [folders, subject, session_num, session_name, locations] = checkSession(folders, subject, session, locations)

[folders, subject] = RCT.dataHandlers.checkSubject(folders, subject);

%% checks session

if(isempty(session) || (isnumeric(session) && (session < 1 || session > numel(folders.sessions.(subject)))))
    error("invalid session.");
end
if(ischar(session))
    subject = string(session);
end

if(isstring(session))
    session_num = find(session == [folders.sessions.(subject)(:).name],1);
    session_name = session;
    if(isempty(session_num))
        error("Session not found.");
    end
else
    session_num = session;
    session_name = folders.sessions.(subject)(session_num).name;
end

%% checks for recording locations
available_locations = folders.sessions.(subject)(session_num).locations;
if(nargin < 4 || isempty(locations))
    locations = available_locations;
elseif(ischar(locations))
    locations = upper(string(locations));
end

if(~all(ismember(locations, available_locations)))
    error("Recording location not found.");
end