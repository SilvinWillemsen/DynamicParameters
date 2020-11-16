stateAtFile = path + mode + "/stateAt" + version + ".csv";
plotIdxFile = path + mode + "/plotIdx" + version + ".csv";
outputFile = path + mode + "/output" + version + ".csv";
cSaveFile = path + mode + "/cSave" + version + ".csv";
NSaveFile = path + mode + "/NSave" + version + ".csv";
NChangeFile = path + mode + "/NChange" + version + ".csv";
lambdaSqSaveFile = path + mode + "/lambdaSqSave" + version + ".csv";
alfTickSaveFile = path + mode + "/alfTickSave" + version + ".csv";
curFsFile = path + mode + "/curFs.csv";

stateAtD = load (stateAtFile);
disp ("stateAt" + version + " loaded")
plotIdxD = load (plotIdxFile);
disp ("plotIdx" + version + " loaded")
outputD = load (outputFile);
disp ("output" + version + " loaded");
cSaveD = load (cSaveFile);
disp ("cSave" + version + " loaded")
NSaveD = load (NSaveFile);
disp ("NSave" + version + " loaded")
NChangeD = load (NChangeFile);
disp ("NChange" + version + " loaded")
lambdaSqSaveD = load (lambdaSqSaveFile);
disp ("LambdaSqSave" + version + " loaded")
alfTickSaveD = load (alfTickSaveFile);
disp ("AlfTickSave" + version + " loaded")