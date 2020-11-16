stateAtFile = path + mode + "/stateAt" + version + ".csv";
plotIdxFile = path + mode + "/plotIdx" + version + ".csv";
outputFile = path + mode + "/output" + version + ".csv";
cSaveFile = path + mode + "/cSave" + version + ".csv";
NSaveFile = path + mode + "/NSave" + version + ".csv";
NChangeFile = path + mode + "/NChange" + version + ".csv";
lambdaSqSaveFile = path + mode + "/lambdaSqSave" + version + ".csv";
alfTickSaveFile = path + mode + "/alfTickSave" + version + ".csv";
curFsFile = path + mode + "/curFs.csv";

stateAtI = load (stateAtFile);
disp ("stateAt" + version + " loaded")
plotIdxI = load (plotIdxFile);
disp ("plotIdx" + version + " loaded")
outputI = load (outputFile);
disp ("output" + version + " loaded");
cSaveI = load (cSaveFile);
disp ("cSave" + version + " loaded")
NSaveI = load (NSaveFile);
disp ("NSave" + version + " loaded")
NChangeI = load (NChangeFile);
disp ("NChange" + version + " loaded")
lambdaSqSaveI = load (lambdaSqSaveFile);
disp ("LambdaSqSave" + version + " loaded")
alfTickSaveI = load (alfTickSaveFile);
disp ("AlfTickSave" + version + " loaded")