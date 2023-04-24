function [Ranks] = getDefaultRanks(R, R2, R3, R4, R5)
if(nargin < 1 || isempty(R))
    R  = 32;
end
if(nargin < 2 || isempty(R2))
    R2 = 32;
end
if(nargin < 3 || isempty(R3))
    R3 = 32;
end
if(nargin < 4 || isempty(R4))
    R4 = 16;
end
if(nargin < 5 || isempty(R5))
    R5 = 8;
end

Ranks.stimulus = R;
Ranks.response = R2;
Ranks.coupling_local = R3;
Ranks.coupling_inter = R4;

Ranks.fixation = R5;

%Ranks.trialwise_latent = R; 
