function [] = setupGMLMpaths()

gmlmPath =  which("kgmlm.GMLM");

if(isempty(gmlmPath))
    addpath ../GMLM/;
    addpath ../GMLM/example;
end