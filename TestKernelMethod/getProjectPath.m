function out = getProjectPath

addpath ../AuxiliaryFunctions/

currentFolder = pwd;

cd ../
out = projectPath;

cd(currentFolder)

end