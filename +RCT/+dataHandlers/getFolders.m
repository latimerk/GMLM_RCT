function [fs] = getFolders()
fs.subjects = "Mingus";

RCT_dir = what("RCT");
fs.RCT_base = fileparts(RCT_dir.path);

RCT_results_base = fs.RCT_base;
RCT_data_base    = fs.RCT_base;

for mm = 1:numel(fs.subjects)
    sub = fs.subjects(mm);

    fs.data.raw.(sub)       = sprintf("%s/Data/raw/%s/", RCT_data_base, sub);
    fs.data.processed.(sub) = sprintf("%s/Data/processed/%s/", RCT_data_base, sub);
    
    fs.results.GMLM.baseline.(sub)      = sprintf("%s/Results/GMLM/baseline/%s/",  RCT_results_base, sub);
    fs.results.GMLM.noCoupling.(sub)    = sprintf("%s/Results/GMLM/Coupling_none/%s/",  RCT_results_base, sub);
    fs.results.GMLM.localCoupling.(sub) = sprintf("%s/Results/GMLM/Coupling_local/%s/", RCT_results_base, sub);
    fs.results.GMLM.allCoupling.(sub)   = sprintf("%s/Results/GMLM/Coupling_all/%s/",   RCT_results_base, sub);
    
    % get session idxs
    raw_fs = dir(fs.data.raw.(sub));
    sessions = [];
    for ii = 1:numel(raw_fs)
        A = split(string(raw_fs(ii).name),["_" "."]);
        for jj = 1:numel(A)
            nn = str2double(A(jj));
            if(~isempty(nn) && ~isnan(nn) && nn ~= 0)
                location = A(jj+1); % why is this only saved in the filename?!
                sessions = cat(1, sessions, [A(jj) location]);

                break;
            end
        end
    end
    if(~isempty(sessions))
        session_names = unique(sessions(:,1));
        fs.sessions.(sub) = struct("name", cell(numel(session_names),1), "locations", []);
        for ii = 1:numel(session_names)
            fs.sessions.(sub)(ii).name = session_names(ii);
            fs.sessions.(sub)(ii).locations = upper(unique(sessions(session_names(ii) == sessions(:,1), 2))');
        end
    else
        warning("No data found for subject %s.", fs.subjects(mm));
    end
end