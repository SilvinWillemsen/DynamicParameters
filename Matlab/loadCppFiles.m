% clear all
% mode = "Release";  % Release or Debug
% version = "D" % (D)ynamic or (I)nterpolation

path = "/Users/SilvinW/Desktop/output/FullGridInterpolation-ejjbtghhbugjpfgqiqqcdpfkqwvo/Build/Products/";

versionFile = path + mode + "/curVersion.txt";
origVersion = fileread(versionFile);
version = origVersion;
if origVersion == "B"
    version = "I";
    loadIFiles;
    version = "D";
    loadDFiles;
elseif origVersion == "I"
     loadIFiles;

elseif origVersion == "D"
     loadDFiles;
end

load (curFsFile)
disp ("curFs loaded")
