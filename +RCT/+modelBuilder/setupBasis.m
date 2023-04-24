%% creates basis functions for modeling the RCT task
% 
%  Uses modified cardinal splines from
%   Sarmashghi, M., Jadhav, S. P., & Eden, U. (2021). Efficient Spline Regression for Neural Spiking Data. bioRxiv.
%
%
function [bases] = setupBasis(binSize_ms, orthogonalize, orthogonalize_spkHist)
if(nargin < 1 || isempty(binSize_ms))
    binSize_ms = 1;
end
if(nargin < 2 || isempty(orthogonalize))
    orthogonalize = true;
end
if(nargin < 3 || isempty(orthogonalize_spkHist))
    orthogonalize_spkHist = true;
end
dd = 1e-12;

%% spike history basis
bases.spkHist.s    = []; % empty = default
if(binSize_ms < 2)
    bases.spkHist.c_pt = [0 2 4 10:10:40 60 80 120 160 220 280 340 400];
elseif(binSize_ms < 4)
    bases.spkHist.c_pt = [0 4 10 20 30 40 60 80 120 160 220 280 340 400];
elseif(binSize_ms < 10)
    bases.spkHist.c_pt = [0 10 20 30 40 60 80 120 160 220 280 340 400];
else
    error("Can't setup spike history for this bin size.");
end
bases.spkHist.c_pt(end) = bases.spkHist.c_pt(end) + dd;
bases.spkHist.lag  = 400;

bases.spkHist.zeroEndPoints = [false false];
bases.spkHist.scale_goal = nan;


[bases.spkHist.B, bases.spkHist.B_0, bases.spkHist.tts_0, bases.spkHist.tts] = DMC.modelBuilder.ModifiedCardinalSpline(bases.spkHist.lag, bases.spkHist.c_pt, bases.spkHist.s, binSize_ms, bases.spkHist.zeroEndPoints);

bases.spkHist.std.center = 1;
bases.spkHist.std.width  = 40;
bases.spkHist.std.decay  = 1/360;
bases.spkHist.std.min    = 0.5;
bases.spkHist.std.mask   = bases.spkHist.tts_0(:) > 5;
bases.spkHist.std.mask_std = 2;
%% stimulus filter basis (same component for noise stim on, motion stim on, targets on, and go signal)

% max time between noise on and saccade is ~1.5 s

bases.stimulus.s    = []; % empty = default
% bases.stimulus.c_pt = [0:10:40 45:5:70 80:10:100 120:20:200 250:50:1600];
bases.stimulus.c_pt = [0:10:40 45:5:120 130:10:200 220:20:300 350:50:1200];
bases.stimulus.window  = [0 1200];
bases.stimulus.c_pt(1  ) = bases.stimulus.c_pt(1) - dd;
bases.stimulus.c_pt(end) = bases.stimulus.c_pt(end) + dd;
bases.stimulus.zeroEndPoints = [true false];
bases.stimulus.scale_goal = 1;


[bases.stimulus.B, bases.stimulus.B_0, bases.stimulus.tts_0, bases.stimulus.tts] = DMC.modelBuilder.ModifiedCardinalSpline(bases.stimulus.window, bases.stimulus.c_pt, bases.stimulus.s, binSize_ms, bases.stimulus.zeroEndPoints);

bases.stimulus.std.center = 100;
bases.stimulus.std.width  = 100;
bases.stimulus.std.decay  = 1/(1200-200);
bases.stimulus.std.min    = 0.25;
bases.stimulus.std.mask   = bases.stimulus.tts_0(:) >= 0;
bases.stimulus.std.mask_std = 1;
%% responses
bases.response.s       = []; % empty = default
% bases.response.c_pt    = [-300 -250 -200 -180 -160:10:-40 -20 0 20 40 60];
% bases.response.c_pt    = [-300:20:-200 -190:10:-140 -135:5:-60 -50:10:60];
bases.response.c_pt    = [-300:20:-200 -190 -180:5:-50 -40:10:60];
bases.response.c_pt(1  ) = bases.response.c_pt(1) - dd;
bases.response.c_pt(end) = bases.response.c_pt(end) + dd;
bases.response.window  = [-300 60];
bases.response.zeroEndPoints = [false false];

[bases.response.B, bases.response.B_0, bases.response.tts_0, bases.response.tts] = DMC.modelBuilder.ModifiedCardinalSpline(bases.response.window, bases.response.c_pt, bases.response.s, binSize_ms, bases.response.zeroEndPoints);

bases.response.std.center = -80;
bases.response.std.width  = 60;
bases.response.std.decay  = 1/160;
bases.response.std.min    = 0.5;
bases.response.std.mask   = bases.response.tts_0(:) <= 60;
bases.response.std.mask_std = 1;

%% fixation (actually aligned to noise_on, but in past. Only used for baseline rate model)
bases.fixation.s       = []; % empty = default
bases.fixation.c_pt    = -400:50:50;
bases.fixation.c_pt(1  ) = bases.fixation.c_pt(1) - dd;
bases.fixation.c_pt(end) = bases.fixation.c_pt(end) + dd;
bases.fixation.window  = [-400 50];
bases.fixation.zeroEndPoints = [false true];

[bases.fixation.B, bases.fixation.B_0, bases.fixation.tts_0, bases.fixation.tts] = DMC.modelBuilder.ModifiedCardinalSpline(bases.fixation.window, bases.fixation.c_pt, bases.fixation.s, binSize_ms, bases.fixation.zeroEndPoints);

bases.fixation.std.center = -50;
bases.fixation.std.width  = 50;
bases.fixation.std.decay  = 1/300;
bases.fixation.std.min    = 1;
bases.fixation.std.mask   = bases.fixation.tts_0(:) <= 40;
bases.fixation.std.mask_std = 1;
%%
fs = fieldnames(bases);

for ii = 1:numel(fs)
    bases.(fs{ii}).tts_idx = 1:size(bases.(fs{ii}).B,1);
    if(~orthogonalize && ~startsWith(fs{ii}, "spk"))
        bases.(fs{ii}).B   = bases.(fs{ii}).B_0;
    elseif(~orthogonalize_spkHist && startsWith(fs{ii}, "spk"))
        bases.(fs{ii}).B   = bases.(fs{ii}).B_0;
    end
end
bases.binSize_ms = binSize_ms;

%%
% _bin_size_ms _c_pt _tension _zero_first _zero_last 