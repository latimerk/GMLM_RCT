fname = "Data/processed/Mingus/Session_20201211_complete_curation.mat";

badCells.("FEF") = [];
badCells.("LIP") = ["20201211_lip_89";
                    "20201211_lip_92";
                    "20201211_lip_110"];
badCells.("SC")  = ["20201211_sc_45";
                    "20201211_sc_169"];

save(fname, "badCells", "-v7.3");