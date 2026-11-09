function cases=qstab_cases_manifestfix(checkManifest)
%QSTAB_CASES_MANIFESTFIX Keep empty path components visible to validation.
% Pure data checks: no filesystem access, solver execution or path changes.
cases={ ...
 'Repeated forward separators are invalid, not just missing source records',@()forwardSeparators(checkManifest); ...
 'Repeated and mixed backslash separators remain invalid after normalization',@()backslashSeparators(checkManifest); ...
 'Empty path components are rejected even in ignored non-MATLAB entries',@()nonSourceEntries(checkManifest); ...
 'A matching malformed inventory cannot turn an invalid manifest into success',@()matchedInventory(checkManifest); ...
 'Empty endpoints and dot-traversal segments remain explicitly invalid',@()endpointsAndDots(checkManifest); ...
 'Canonical slash variants, dotted names, comments and line endings remain valid',@()canonicalPaths(checkManifest)};
end

function [text,I]=fixture()
h=repmat('a',1,64);g=repmat('b',1,64);
text=sprintf('%s  src/+qstab/eig.m\n%s  examples/demo.m\n',h,g);
s=struct('path','src/+qstab/eig.m','sha256',h,'bytes',3,'error','');
s(2)=struct('path','examples/demo.m','sha256',g,'bytes',4,'error','');
I=struct('complete',true,'sources',s);
end

function rejectPaths(f,paths)
[t,I]=fixture();saved=I;
for ignoreCase=[false true]
    for k=1:numel(paths)
        line=sprintf('%s  %s\n',repmat('c',1,64),paths{k});
        R=f([t line],I,ignoreCase);
        assert(~R.ok&&~R.manifestValid, ...
            'qstab:ManifestFix:AcceptedPath','Unsafe manifest path accepted: %s',paths{k});
        assert(R.expectedCount==2&&numel(R.errors)==1);
        assert(~isempty(strfind(R.errors{1},'Unsafe manifest path'))); %#ok<STREMP>
        assert(isempty(R.missing)&&isempty(R.modified)&&isempty(R.unverified));
        assert(isequaln(I,saved));
    end
end
end

function forwardSeparators(f)
rejectPaths(f,{'src//evil.m','src///evil.m','src/+qstab//eig.m', ...
    'src//+qstab/eig.m','examples////demo.m'});
end

function backslashSeparators(f)
b=char(92);
rejectPaths(f,{['src' b b 'evil.m'],['src' b b b 'evil.m'], ...
    ['src/' b 'evil.m'],['src' b '/evil.m'], ...
    ['src' b '+qstab' b b 'eig.m']});
end

function nonSourceEntries(f)
b=char(92);
rejectPaths(f,{'docs//README.md','docs///data.json','verification//data.txt', ...
    'docs//LICENSE',['docs' b b 'README.md'],['docs/' b 'README.md']});
end

function matchedInventory(f)
[t,I]=fixture();bad='src//evil.m';h=repmat('c',1,64);
I.sources(3)=struct('path',bad,'sha256',h,'bytes',7,'error','');
for ignoreCase=[false true]
    R=f([t sprintf('%s  %s\n',h,bad)],I,ignoreCase);
    assert(~R.ok&&~R.manifestValid&&R.expectedCount==2);
    assert(isequal(R.unexpected,{bad}));
    assert(numel(R.errors)==1&&~isempty(strfind(R.errors{1},'Unsafe manifest path'))); %#ok<STREMP>
end
end

function endpointsAndDots(f)
rejectPaths(f,{'/evil.m','//evil.m','src/','src//','src/./evil.m', ...
    'src/../evil.m','./evil.m','../evil.m','src/..//evil.m'});
end

function canonicalPaths(f)
[t,I]=fixture();h=repmat('c',1,64);
good={'src/+qstab/+internal/core.m','examples/review case.m','src/name..backup.m'};
for k=1:numel(good)
    t=[t sprintf('%s  %s\n',h,good{k})]; %#ok<AGROW>
    I.sources(end+1)=struct('path',good{k},'sha256',h,'bytes',k,'error',''); %#ok<AGROW>
end
for ignoreCase=[false true]
    R=f(t,I,ignoreCase);assert(R.ok&&R.manifestValid&&R.expectedCount==5);
    windowsText=strrep(t,'/',char(92));
    windowsText=strrep(windowsText,sprintf('\n'),sprintf('\r\n'));
    R=f([sprintf('# comment\r\n\r\n') windowsText],I,ignoreCase);
    assert(R.ok&&R.manifestValid&&R.expectedCount==5);
end
end
