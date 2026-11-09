function cases=qstab_cases_namespace(root)
%QSTAB_CASES_NAMESPACE Public/private separation without invoking user aliases.
cases={ ...
    'Public distribution contains no q-prefixed compatibility aliases',@noAliases; ...
    'A foreign qeig stays separate from the public right-spectrum namespace',@separateEigen; ...
    'A foreign qisstable is never called by the public nominal solver',@foreignStable; ...
    'Removed wrapper opt-in is rejected before any path changes',@rejectOptIn; ...
    'Public setup preserves foreign qfiles paths and contents',@preserveForeign};
    function noAliases()
        assert(exist(fullfile(root,'qwrappers'),'dir')~=7);
        assert(exist(fullfile(root,'qeig.m'),'file')~=2);
        assert(exist(fullfile(root,'src','+qstab','eig.m'),'file')==2);
    end
    function separateEigen()
        [folder,guard]=fixture('qeig'); %#ok<ASGLU,NASGU>
        addpath(folder,'-begin');resolved=which('qeig');
        [z,I]=qstab.eig(diag([-1,-2]));
        assert(I.isValid&&isequal(sort(real(z)),[-2;-1]));
        assert(strcmp(which('qeig'),resolved));
    end
    function foreignStable()
        [folder,guard]=fixture('qisstable'); %#ok<ASGLU,NASGU>
        addpath(folder,'-begin');before=path;
        [~,R] = qstab.isstable([2 3 1]);assert(strcmp(R.status,'stable'));
        assert(strcmp(path,before));
    end
    function rejectOptIn()
        before=path;rejected=false;
        try,qstab_setup('Wrappers',true);catch ME,rejected=strcmp(ME.identifier,'qstab:setup:PrivateWrappers');end
        assert(rejected&&strcmp(path,before));
    end
    function preserveForeign()
        [folder,guard]=fixture('qisstable'); %#ok<ASGLU,NASGU>
        addpath(folder,'-end');bytes=fileread(fullfile(folder,'qisstable.m'));
        I=qstab_setup('Wrappers',false,'Examples',false);
        assert(~I.wrappersEnabled&&any(strcmp(strsplit(path,pathsep),folder)));
        assert(strcmp(bytes,fileread(fullfile(folder,'qisstable.m'))));
    end
end
function [folder,guard]=fixture(name)
original=path;folder=tempname;mkdir(folder);
guard=onCleanup(@()cleanup(original,folder));
fid=fopen(fullfile(folder,[name '.m']),'wt');assert(fid>=0);
fprintf(fid,'function varargout=%s(varargin)\nerror(''qstab:ForeignExecuted'',''Foreign alias must never run.'');\nend\n',name);
fclose(fid);rehash;
end
function cleanup(original,folder)
path(original);if exist(folder,'dir')==7,rmdir(folder,'s');end
end
