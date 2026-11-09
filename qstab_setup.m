function info=qstab_setup(varargin)
%QSTAB_SETUP Connect the public namespace and optional examples for this session.
%   INFO=QSTAB_SETUP() adds this package root, src, and examples to path.
%   INFO=QSTAB_SETUP('Examples',false) does not add the examples directory.
%   QSTAB_SETUP without an output also prints the version and fixed report ZIP.
%
%   Public qstab contains no global q-prefixed wrappers. These belong only to
%   the separate private qfiles integration. For compatibility, Wrappers=false
%   is accepted; Wrappers=true raises qstab:setup:PrivateWrappers BEFORE any
%   path changes. No private toolbox is loaded or initialized.
%
%   Only the parent src of +qstab is added. Never add +qstab directly, and do
%   not use genpath on the whole distribution. This function does not modify
%   startup.m, savepath, global variables, existing qfiles, or class definitions.
%   A stale qwrappers path inside THIS installation is disconnected but its
%   files are not deleted; the installation checker reports leftover sources.
%
%   Example:
%       % Start in the extracted qstab_public directory.
%       info=qstab_setup('Examples',true); % silent setup
%       [~,R] = qstab.isstable([2 3 1]);
%       disp(R.status);
%
%   Result archive: <this root>/qstab_test_results.zip, independent of pwd.
%   See also run_qstab_public, test_qstab, qstab_check_installation, qstab.about.
ip=inputParser;
addParameter(ip,'Wrappers',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Examples',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});
if ip.Results.Wrappers
    error('qstab:setup:PrivateWrappers', ...
        'q-prefixed aliases are private qfiles adapters and are not distributed in public qstab. Use namespace calls. No path change was made.');
end
root=fileparts(mfilename('fullpath'));
if exist(fullfile(root,'src','+qstab','isstable.m'),'file')~=2|| ...
        (ip.Results.Examples&&exist(fullfile(root,'examples'),'dir')~=7)
    error('qstab:setup:IncompleteInstallation','Required source/example directory is missing. Re-extract the complete distribution.');
end
addpath(root);addpath(fullfile(root,'src'));
if ip.Results.Examples,addpath(fullfile(root,'examples'));end
wrapperDir=fullfile(root,'qwrappers');dirs=strsplit(path,pathsep);
if ispc,present=any(strcmpi(dirs,wrapperDir));else,present=any(strcmp(dirs,wrapperDir));end
if present,rmpath(wrapperDir);end
info=struct('root',root,'namespaceParent',fullfile(root,'src'), ...
    'wrappersEnabled',false,'examplesEnabled',ip.Results.Examples, ...
    'version','0.18.0','wrapperPolicy','private adapters are distributed separately; public namespace only');
if nargout==0
    fprintf('qstab 0.18.0 ready: %s\n',root);
    fprintf('Next: report = run_qstab_public;\n');
    fprintf('Result ZIP (always): %s\n',fullfile(root,'qstab_test_results.zip'));
end
end
