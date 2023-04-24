function [S] = getSTDForKernels(tts, settings)


S = (1-settings.min)*max(0, 1 - settings.decay*max(0, abs(tts - settings.center) - settings.width))  + settings.min;