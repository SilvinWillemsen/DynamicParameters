stateAtFile = path + mode + "/stateAt" + version + ".csv";
plotIdxFile = path + mode + "/plotIdx" + version + ".csv";
outputFile = path + mode + "/output" + version + ".csv";
cSaveFile = path + mode + "/cSave" + version + ".csv";
NSaveFile = path + mode + "/NSave" + version + ".csv";
NChangeFile = path + mode + "/NChange" + version + ".csv";
lambdaSqSaveFile = path + mode + "/lambdaSqSave" + version + ".csv";
alfTickSaveFile = path + mode + "/alfTickSave" + version + ".csv";
aUSaveFile = path + mode + "/aUSave" + version + ".csv";
curFsFile = path + mode + "/curFs.csv";
diffSaveFile = path + mode + "/diffSave.csv";
timeIdxFile = path + mode + "/timeIdxSave.csv";


if ~onlyOutput
    stateAtD = load (stateAtFile);
    disp ("stateAt" + version + " loaded")
    plotIdxD = load (plotIdxFile);
    disp ("plotIdx" + version + " loaded")
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
    aUSaveD = load (aUSaveFile);
    disp ("AUSave" + version + " loaded")
    diffSaveD = load (diffSaveFile);
    disp ("diffSave" + version + " loaded")
    timeIdxD = load (timeIdxFile);
    disp ("timeIdx" + version + " loaded")

end
outputD = load (outputFile);
disp ("output" + version + " loaded");
