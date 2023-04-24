
% optional key/value pairs are:
%  includeDirection  = true/false (DEFAULT: true) cosine direction tuning only for now
%  includeCategory   = true/false (DEFAULT: true)
%  includeTargetLocationsOnset  = true/false (DEFAULT: true)
%  includeTargetLocationsFixOffset = true/false (DEFAULT: false)
%  includeResponseLocations   = true/false (DEFAULT: true)
%  includeResponseTargetColor = true/false (DEFAULT: true)
%  includeResponseLocationTargetInteractions = true/false (DEFAULT: true)
%  includeFixation  = true/false (DEFAULT: false) Include event for fixation time

function [stimConfig, responseConfig, fixationConfig] = getTaskModelSetup(TaskInfo, varargin)

p = inputParser;
p.CaseSensitive = false;

addRequired( p, "TaskInfo"     ,     @isstruct);
addParameter(p, "includeDirection"    ,  true,    @islogical);
addParameter(p, "includeCategory"  ,  true,    @islogical);
addParameter(p, "includeTargetLocationsOnset"  ,  true,    @islogical);
addParameter(p, "includeTargetLocationsFixOffset" ,  true,   @islogical);
addParameter(p, "includeResponseLocations"  ,  true,    @islogical);
addParameter(p, "includeResponseTargetColor"  ,  true,    @islogical);
addParameter(p, "includeResponseLocationTargetInteractions"  ,  true,    @islogical);
addParameter(p, 'correctSTDOfGroups',  true,    @islogical);
addParameter(p, 'coefficientSTD', 1, @isnumeric);



parse(p, TaskInfo, varargin{:});
% then set/get all the inputs out of this structure
includeDirection  = p.Results.includeDirection;
includeCategory   = p.Results.includeCategory;
includeTargetLocationsOnset  = p.Results.includeTargetLocationsOnset;
includeTargetLocationsFixOffset = p.Results.includeTargetLocationsFixOffset;
includeResponseLocations   = p.Results.includeResponseLocations;
includeResponseTargetColor = p.Results.includeResponseTargetColor;
includeResponseLocationTargetInteractions = p.Results.includeResponseLocationTargetInteractions;
coefficientSTD = p.Results.coefficientSTD;
correctSTDOfGroups = p.Results.correctSTDOfGroups;

%%
unique_cats = unique(TaskInfo.categories);
NC = numel(unique_cats);
if(NC ~= 2)
    error("This function only handles 2 category tasks.");
end
ND = numel(TaskInfo.directions);

%% setup basic stim info

% events: noise stim, dots stim, targets on, targets off (go)

totalParams = 4; % event onsets
totalParams = totalParams + includeTargetLocationsOnset; % green/red location for targets on
totalParams = totalParams + includeTargetLocationsFixOffset; % green/red location for targets off
totalParams = totalParams + includeCategory; % category of noise direction
totalParams = totalParams + includeDirection*2; % cosine/sine of noise direction

stimConfig = struct("event_name", cell(4,1), "coeff_names", [], "coeff_groups", [], "regressors", []);

ctr = 0;

ss = 1;
stimConfig(ss).event_name = "noise_on";
stimConfig(ss).coeff_names = "noise_on";
stimConfig(ss).coeff_groups = {1 + ctr};
ctr = ctr + 1;
stimConfig(ss).regressors = zeros(1, totalParams);
stimConfig(ss).regressors(1, stimConfig(ss).coeff_groups{1}) = 1;

ss = ss + 1;
stimConfig(ss).event_name = "stim_on";
stimConfig(ss).coeff_names = "dots_on";
stimConfig(ss).coeff_groups = {1 + ctr};
ctr = ctr + 1;
stimConfig(ss).regressors = zeros(ND, totalParams);
stimConfig(ss).regressors(:, stimConfig(ss).coeff_groups{1}) = 1;

if(includeCategory)
    stimConfig(ss).coeff_names  = [stimConfig(ss).coeff_names "dots_cat"];
    stimConfig(ss).coeff_groups = cat(2,stimConfig(ss).coeff_groups, {1 + ctr});
    stimConfig(ss).regressors(TaskInfo.categories == unique_cats(1), stimConfig(ss).coeff_groups{end}) =  1;
    stimConfig(ss).regressors(TaskInfo.categories == unique_cats(2), stimConfig(ss).coeff_groups{end}) = -1;
    ctr = ctr + 1;
end
if(includeCategory)
    stimConfig(ss).coeff_names  = [stimConfig(ss).coeff_names "dots_sin_cos_dir"];
    stimConfig(ss).coeff_groups = cat(2,stimConfig(ss).coeff_groups, {[1 2] + ctr});
    stimConfig(ss).regressors(:, stimConfig(ss).coeff_groups{end}(1)) = sind(TaskInfo.directions);
    stimConfig(ss).regressors(:, stimConfig(ss).coeff_groups{end}(2)) = cosd(TaskInfo.directions);
    ctr = ctr + 2;
end

ss = ss + 1;
stimConfig(ss).event_name = "targets_on";
stimConfig(ss).coeff_names = "targets_on";
stimConfig(ss).coeff_groups = {1 + ctr};
ctr = ctr + 1;
stimConfig(ss).regressors = zeros(2, totalParams);
stimConfig(ss).regressors(:, stimConfig(ss).coeff_groups{1}) = 1;
if(includeTargetLocationsOnset)
    stimConfig(ss).coeff_names  = [stimConfig(ss).coeff_names "targets_rg_locations"];
    stimConfig(ss).coeff_groups = cat(2,stimConfig(ss).coeff_groups, {1 + ctr});
    stimConfig(ss).regressors(1, stimConfig(ss).coeff_groups{end}) =  1;
    stimConfig(ss).regressors(2, stimConfig(ss).coeff_groups{end}) = -1;
    ctr = ctr + 1;
end

ss = ss + 1;
stimConfig(ss).event_name = "fix_off";
stimConfig(ss).coeff_names = "fix_off";
stimConfig(ss).coeff_groups = {1 + ctr};
ctr = ctr + 1;
stimConfig(ss).regressors = zeros(2, totalParams);
stimConfig(ss).regressors(:, stimConfig(ss).coeff_groups{1}) = 1;
if(includeTargetLocationsFixOffset)
    stimConfig(ss).coeff_names  = [stimConfig(ss).coeff_names "fix_off_rg_locations"];
    stimConfig(ss).coeff_groups = cat(2,stimConfig(ss).coeff_groups, {1 + ctr});
    stimConfig(ss).regressors(1, stimConfig(ss).coeff_groups{end}) =  1;
    stimConfig(ss).regressors(2, stimConfig(ss).coeff_groups{end}) = -1;
    ctr = ctr + 1;
end


%%
totalParams_r = 1 + includeResponseLocations + includeResponseTargetColor + includeResponseLocationTargetInteractions*2;
responseConfig = struct("event_name", cell(1,1), "coeff_names", [], "coeff_groups", [], "regressors", []);
ctr = 0;

ss = 1;
responseConfig(ss).event_name = "saccade_end";
responseConfig(ss).coeff_names = "saccade";
responseConfig(ss).coeff_groups = {1 + ctr};
ctr = ctr + 1;
responseConfig(ss).regressors = zeros(4, totalParams_r);
responseConfig(ss).regressors(:, responseConfig(ss).coeff_groups{1}) = 1;
    %rows: 1  = left/down  red
    %      2  = left/down  green
    %      3  = right/up red
    %      4  = right/up green

if(includeResponseLocations)
    responseConfig(ss).coeff_names  = [responseConfig(ss).coeff_names "saccade_location"];
    responseConfig(ss).coeff_groups = cat(2,responseConfig(ss).coeff_groups, {1 + ctr});
    responseConfig(ss).regressors(1:2, responseConfig(ss).coeff_groups{end}) =  1;
    responseConfig(ss).regressors(3:4, responseConfig(ss).coeff_groups{end}) = -1;
    ctr = ctr + 1;
end
if(includeResponseTargetColor)
    responseConfig(ss).coeff_names  = [responseConfig(ss).coeff_names "saccade_color"];
    responseConfig(ss).coeff_groups = cat(2,responseConfig(ss).coeff_groups, {1 + ctr});
    responseConfig(ss).regressors([1 3], responseConfig(ss).coeff_groups{end}) =  1;
    responseConfig(ss).regressors([2 4], responseConfig(ss).coeff_groups{end}) = -1;
    ctr = ctr + 1;
end
if(includeResponseLocationTargetInteractions)
    responseConfig(ss).coeff_names  = [responseConfig(ss).coeff_names "saccade_color_loc1" "saccade_color_loc2"];
    responseConfig(ss).coeff_groups = cat(2,responseConfig(ss).coeff_groups, {1 + ctr}, {2 + ctr});
    responseConfig(ss).regressors(1, responseConfig(ss).coeff_groups{end}-1) =   1;
    responseConfig(ss).regressors(2, responseConfig(ss).coeff_groups{end}-1) =  -1;
    responseConfig(ss).regressors(3, responseConfig(ss).coeff_groups{end}  ) =   1;
    responseConfig(ss).regressors(4, responseConfig(ss).coeff_groups{end}  ) =  -1;
    ctr = ctr + 2;
end

%% setup fixation info: aligned to noise on

% events: noise stim
totalParams = 1;

fixationConfig = struct("event_name", cell(1,1), "coeff_names", [], "coeff_groups", [], "regressors", []);

ctr = 0;

ss = 1;
fixationConfig(ss).event_name = "noise_on";
fixationConfig(ss).coeff_names = "noise_on";
fixationConfig(ss).coeff_groups = {1 + ctr};
ctr = ctr + 1;
fixationConfig(ss).regressors = zeros(1, totalParams);
fixationConfig(ss).regressors(1, fixationConfig(ss).coeff_groups{1}) = 1;

%% NORMALIZATION
% like z-scoring. Makes the shrinkage comparable across the different coefficients.
if(correctSTDOfGroups)
    for ss = 1:numel(stimConfig)
        stimConfig(ss).coeff_group_normalizer = nan(numel(stimConfig(ss).coeff_groups), 1);
        for gg = 1:numel(stimConfig(ss).coeff_groups)
            gg_idx = stimConfig(ss).coeff_groups{gg};
            if(~isempty(gg_idx))
                %sig = std(stimConfig(ss).regressors(:, gg_idx), 1, "all");
                sig2 = mean(stimConfig(ss).regressors(:, gg_idx).^2, 1); %% assumes 0 mean
                if(all(sig2 > 0))
                    cc = coefficientSTD./sqrt(sum(sig2));
                    stimConfig(ss).coeff_group_normalizer(gg) = cc;
                    stimConfig(ss).regressors(:, gg_idx) = stimConfig(ss).regressors(:, gg_idx)*cc;
                end
            end
        end
    end
    for ss = 1:numel(responseConfig)
        responseConfig(ss).coeff_group_normalizer = nan(numel(responseConfig(ss).coeff_groups), 1);
        for gg = 1:numel(responseConfig(ss).coeff_groups)
            gg_idx = responseConfig(ss).coeff_groups{gg};
            if(~isempty(gg_idx))
                %sig = std(responseConfig(ss).regressors(:, gg_idx), 1, "all");
                sig2 = mean(responseConfig(ss).regressors(:, gg_idx).^2, 1); %% assumes 0 mean
                if(all(sig2 > 0))
                    cc = coefficientSTD./sqrt(sum(sig2));
                    responseConfig(ss).coeff_group_normalizer(gg) = cc;
                    responseConfig(ss).regressors(:, gg_idx) = responseConfig(ss).regressors(:, gg_idx)*cc;
                end
            end
        end
    end
    for ss = 1:numel(fixationConfig)
        fixationConfig(ss).coeff_group_normalizer = nan(numel(fixationConfig(ss).coeff_groups), 1);
        for gg = 1:numel(fixationConfig(ss).coeff_groups)
            gg_idx = fixationConfig(ss).coeff_groups{gg};
            if(~isempty(gg_idx))
                %sig = std(responseConfig(ss).regressors(:, gg_idx), 1, "all");
                sig2 = mean(fixationConfig(ss).regressors(:, gg_idx).^2, 1); %% assumes 0 mean
                if(all(sig2 > 0))
                    cc = coefficientSTD./sqrt(sum(sig2));
                    fixationConfig(ss).coeff_group_normalizer(gg) = cc;
                    fixationConfig(ss).regressors(:, gg_idx) = fixationConfig(ss).regressors(:, gg_idx)*cc;
                end
            end
        end
    end
end