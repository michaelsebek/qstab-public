function [value,S]=qstab_public_session(root,callback)
%QSTAB_PUBLIC_SESSION Run a callback with only this package and matlabroot.
% Restores the original path and current folder on success or exception.
% Callback handles are used only in memory and never placed in a saved report.
if ~(isa(callback,'function_handle')&&ischar(root)&&isrow(root))
    error('qstab:PublicSessionInput','Expected a package path and callback.');
end
[found,f]=fileattrib(root);
if ~found||~isscalar(f)||~f.directory
    error('qstab:PublicSessionInput','Package folder does not exist.');
end
root=f.Name;oldPath=path;oldFolder=pwd;
S=qstab_public_path_plan(oldPath,root,matlabroot);
S.originalFolder=oldFolder;
restore=onCleanup(@()localRestore(oldPath,oldFolder));
cd(root);
path(strjoin(S.isolatedPaths,pathsep));
S.activeFolder=pwd;S.activePath=path;
value=callback();
clear restore
S.pathRestored=strcmp(path,oldPath);S.folderRestored=strcmp(pwd,oldFolder);
end
function localRestore(oldPath,oldFolder)
% Functions from matlabroot remain available throughout isolation/cleanup.
path(oldPath);
cd(oldFolder);
end
