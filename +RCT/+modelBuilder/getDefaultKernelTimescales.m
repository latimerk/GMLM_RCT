function [kernelTimescales] = getDefaultKernelTimescales()
L = 20;
L2 = 10;
kernelTimescales.stimulus = L; %IN ms
kernelTimescales.response = L;
kernelTimescales.spkHist = L2;

kernelTimescales.fixation = 50;
end

