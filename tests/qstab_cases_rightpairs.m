function T=qstab_cases_rightpairs(writer)
%QSTAB_CASES_RIGHTPAIRS Independent class counts AND quaternionic eigenbases.
C=qstab_eigen_cases();T=cell(0,2);
for k=1:numel(C)
    id=C(k).id;
    T(end+1,:)={['Right eigenpairs preserve the analytic fixture: ' id],@()fixture(id)}; %#ok<AGROW>
end
T=[T;{ ...
 'Repeated real eigenvectors span all quaternionic coordinates',@realBasis; ...
 'Canonicalizing negative imaginary representatives also rotates the vectors',@negativeAxes; ...
 'Small nonzero imaginary parts are retained in the eigenvalue list',@tinyImaginary; ...
 'Close real eigenvalues are not replaced by a unique representative',@closeReal; ...
 'Dense quaternion similarity preserves repeated real eigenvectors',@()quaternionSimilarity(false); ...
 'Dense quaternion similarity preserves repeated nonreal eigenvectors',@()quaternionSimilarity(true); ...
 'Matrix and eigenvector common-scale checks cover extreme finite inputs',@scales; ...
 'Numeric complex and explicit-component inputs have the same semantics',@numericInput; ...
 'The established eig second output remains diagnostic rather than eigenvectors',@oldContract; ...
 'A single eigpairs output is V and never the eigenvalue vector',@singleOutput; ...
 'Empty input retains zero counts and zero-sized component outputs',@emptyInput; ...
 'Nonsquare and genuinely multidimensional eigenpair inputs are rejected',@badShape; ...
 'Nonfinite quaternion components cannot reach the eigensolver',@nonfinite; ...
 'Invalid output requests fail rather than selecting an arbitrary convention',@badOptions; ...
 'Common scaling cannot silently erase a nonzero input component',@scalingUnderflow; ...
 'Pair diagnostics use each quaternion column and a separate basis rank',@diagnostics; ...
 'Right-eigenpair calls preserve input, path, directory and random state',@unchangedState; ...
 'Standalone right eigenvectors do not call a foreign qeig',@noLegacyDependency; ...
 'Optional primary-vector comparison checks all fourteen references',@comparisonAll; ...
 'The default comparison remains spectrum-only unless vectors are requested',@comparisonDefault; ...
 'Defective primary eigenvectors are retained without claiming a full basis',@comparisonDefective; ...
 'An explicit synthetic legacy solver is kept distinct from primary eigenpairs',@comparisonSynthetic; ...
 'Unknown primary-vector comparison options are rejected',@comparisonOptions; ...
 'Right-eigenpair data and optional comparison checks survive the fixed ZIP',@()roundtrip(writer)}];
end
function fixture(id)
C=qstab_eigen_cases();F=C(strcmp({C.id},id));
check(F.matrix,F.reference,F.diagonalizable);
end
function [V,D,I]=check(A,ref,diagonalizable)
[V,D,I]=qstab.eigpairs(A,'Output','components');
S=qstab_eigen_assess(A,ref,I.standardEigenvalues,'CheckVectors',true, ...
    'Vectors',V,'Diagonal',D,'Info',I,'Diagonalizable',logical(diagonalizable));
assert(S.ok,S.message);assert(I.eigenpairsValid&&I.spectrumValid);
assert(numel(I.standardEigenvalues)==numel(ref)&&I.algebraicCount==numel(ref));
assert(all(abs(I.columnNorms-1)<1e-12));
if diagonalizable,assert(I.fullBasisObserved);end
assert(all(imag(I.standardEigenvalues)>=0));
end
function realBasis()
[V,D,I]=check(-eye(4),-ones(4,1),true); %#ok<ASGLU>
assert(I.quaternionBasisRank==4&&I.complexBasisRank==8);
assert(norm(qstab.chi(-eye(4))*qstab.chi(V)-qstab.chi(V)*qstab.chi(D),'fro')==0);
end
function negativeAxes()
A=comp(-eye(3),diag([-1 -2 -3]),zeros(3),zeros(3));
[~,D,I]=check(A,-1+1i*(1:3).',true);
assert(isequal(diag(D.b),(1:3).')&&I.eigenpairsValid);
end
function tinyImaginary()
A=comp(diag([-1 -2]),diag([1e-14 2e-14]),zeros(2),zeros(2));
[~,~,I]=check(A,[-1+1e-14i;-2+2e-14i],true);
assert(all(imag(I.standardEigenvalues)>0));
assert(max(abs(sort(imag(I.standardEigenvalues))-[1e-14;2e-14]))<1e-27);
end
function closeReal()
e=[-1;-1+1e-11;-1+2e-11];[~,D,I]=check(diag(e),e,true);
assert(numel(unique(diag(D.a)))==3&&I.algebraicCount==3);
end
function quaternionSimilarity(nonreal)
n=4;Z=zeros(n);
S=comp(eye(n)+.1*reshape(1:16,n,n), ...
 .05*[0 1 -1 0;1 0 0 -1;0 1 0 1;1 0 1 0], ...
 .03*[1 0 1 0;0 -1 0 1;1 0 -1 0;0 1 0 1], ...
 .02*[0 1 0 1;-1 0 1 0;0 -1 0 1;-1 0 -1 0]);
if nonreal
 D=comp(diag([-1 -1 -1 -2]),diag([1 0 0 0]),diag([0 1 0 0]),diag([0 0 1 0]));
 ref=[-1+1i;-1+1i;-1+1i;-2];
else
 D=comp(diag([-1 -1 -2 -3]),diag([0 0 1 0]),diag([0 0 0 2]),Z);
 ref=[-1;-1;-2+1i;-3+2i];
end
X=qstab.chi(S)*qstab.chi(D)/qstab.chi(S);
A=comp(real(X(1:n,1:n)),imag(X(1:n,1:n)),real(X(1:n,n+1:end)),imag(X(1:n,n+1:end)));
check(A,ref,true);
end
function scales()
C=qstab_eigen_cases();F=C(strcmp({C.id},'triangular'));
for s=[1e-200,1,1e200]
 A=F.matrix;f=fieldnames(A);for k=1:4,A.(f{k})=s*A.(f{k});end
 [~,~,I]=check(A,s*F.reference,true);assert(max(I.columnResiduals)<1e-12);
end
end
function numericInput()
A=[1+2i .25-1i;0 -2-3i];Z=zeros(2);
[V,D,I]=qstab.eigpairs(A);
[W,E,J]=qstab.eigpairs(comp(real(A),imag(A),Z,Z));
assert(isequaln(V,W)&&isequaln(D,E)&&isequaln(I,J));
check(A,[1+2i;-2+3i],true);
end
function oldContract()
[e,I]=qstab.eig(-eye(3));[V,D,J]=qstab.eigpairs(-eye(3)); %#ok<ASGLU>
assert(isnumeric(e)&&isequal(size(e),[3 1])&&isstruct(I)&&isfield(I,'fullSpectrum'));
assert(isstruct(V)&&isfield(V,'a')&&isequal(size(V.a),[3 3]));
assert(isequal(e,J.standardEigenvalues));
end
function singleOutput()
A=-eye(3);V=qstab.eigpairs(A);[W,~,~]=qstab.eigpairs(A);
assert(isequaln(V,W)&&isequal(size(V.a),[3 3]));
end
function emptyInput()
[V,D,I]=qstab.eigpairs(zeros(0));
assert(isequal(size(V.a),[0 0])&&isequal(size(D.a),[0 0]));
assert(I.isValid&&I.fullBasisObserved&&I.quaternionOrder==0&&I.complexBasisRank==0);
assert(isequal(size(I.standardEigenvalues),[0 1]));
end
function badShape()
throws(@()qstab.eigpairs(ones(2,3)));throws(@()qstab.eigpairs(ones(2,2,2)));
throws(@()qstab.eigpairs([1 2 3]));
end
function nonfinite()
throws(@()qstab.eigpairs(NaN));throws(@()qstab.eigpairs(Inf));
throws(@()qstab.eigpairs(comp(1,0,NaN,0)));
end
function badOptions()
throws(@()qstab.eigpairs(1,'Output','unknown'));
throws(@()qstab.eigpairs(1,'Output',[]));
throws(@()qstab.eigpairs(1,'Unexpected',true));
end
function scalingUnderflow()
try
 qstab.eigpairs(diag([1e300 1e-300]));
 error('qstab:test:ExpectedError','An erased input component must be reported.');
catch ME
 assert(strcmp(ME.identifier,'rightpairs:ScalingUnderflow'));
end
end
function diagnostics()
[V,D,I]=qstab.eigpairs(-eye(3)+diag([1 1],1)); %#ok<ASGLU>
assert(I.eigenpairsValid&&I.spectrumValid&&~I.fullBasisObserved);
assert(strcmp(I.status,'rankDeficient'));
assert(numel(I.columnResiduals)==3&&max(I.columnResiduals)<1e-12);
assert(I.complexBasisRank<6&&I.basisReciprocalCondition<1e-12);
end
function unchangedState()
C=qstab_eigen_cases();A=C(1).matrix;before=A;p=path;d=pwd;r=rng;
qstab.eigpairs(A);assert(isequaln(A,before)&&strcmp(path,p)&&strcmp(pwd,d)&&isequaln(rng,r));
end
function noLegacyDependency()
folder=tempname;mkdir(folder);oldPath=path;
tidy=onCleanup(@()cleanup(folder,oldPath)); %#ok<NASGU>
fid=fopen(fullfile(folder,'qeig.m'),'w');assert(fid>=0);
fprintf(fid,'function varargout=qeig(varargin)\nerror(''qstab:fixture:ForbiddenLegacy'',''Must not call qeig.'');\nend\n');fclose(fid);
addpath(folder,'-begin');check(-eye(2),[-1;-1],true);
end
function comparisonAll()
E=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false,'CheckPrimaryVectors',true);
assert(E.ok&&E.qstabOK&&E.primaryEigenpairsRequested&&E.primaryEigenpairsOK&&E.matchedPrimary==14);
for k=1:14,assert(E.cases(k).primary.eigenpairs.ok);end
end
function comparisonDefault()
E=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false,'Cases','equal_real');
assert(E.ok&&~E.primaryEigenpairsRequested&&~isfield(E.cases.primary,'eigenpairs'));
end
function comparisonDefective()
E=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false,'Cases','jordan','CheckPrimaryVectors',true);
S=E.cases.primary.eigenpairs;assert(E.ok&&S.ok&&~S.vectors.fullBasisObserved&&~S.vectors.fullBasisRequired);
end
function comparisonSynthetic()
E=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false,'Cases', ...
 {'equal_real','repeated_real','mixed_multiplicity'},'CheckPrimaryVectors',true, ...
 'LegacyFunction',@syntheticCorrect);
assert(E.ok&&E.qstabOK&&E.primaryEigenpairsOK&&E.legacyOK&&E.allClassMultisetsAgree);
end
function [V,D]=syntheticCorrect(A)
% Synthetic interface fixture, NOT the installed private qeig.
[V,D,I]=qstab.eigpairs(A);if nargout<=1,V=I.standardEigenvalues;end
end
function comparisonOptions()
throws(@()run_qstab_eigen_comparison('CheckPrimaryVectors','true','Verbose',false));
end
function roundtrip(writer)
base=tempname;mkdir(base);tidy=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
c=struct('group','fixture','name','right-eigenpair roundtrip','status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
 'environment',struct('version','fixture','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','20000101T000000','runId','rightpairs-fixture','nativeMATLABRun',false, ...
 'nativeQuaternionOK',false,'ok',true,'status','OK','passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',c);
r.eigenComparison=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false, ...
 'Cases',{'equal_real','repeated_real'},'CheckPrimaryVectors',true);
r=writer(r,base,'qstab_test_results');dest=fullfile(base,'unpacked');mkdir(dest);unzip(r.archiveFile,dest);
s=load(fullfile(dest,'report.mat'),'report');assert(isequaln(s.report.eigenComparison,r.eigenComparison));
f=dir(dest);f=f(~[f.isdir]);assert(isequal(sort({f.name}),sort({'report.mat','report.txt'})));
end
function cleanup(folder,oldPath)
path(oldPath);if exist(folder,'dir')==7,rmdir(folder,'s');end
end
function throws(f)
failed=false;try,f();catch,failed=true;end;assert(failed);
end
function C=comp(a,b,c,d)
C=struct('a',a,'b',b,'c',c,'d',d);
end
