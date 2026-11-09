function T=qstab_cases_eigencomparison(writer)
%QSTAB_CASES_EIGENCOMPARISON Diagnostic regressions, not a replacement qeig.
T={ ...
 'Eigen comparison catalogue has fourteen deterministic analytic fixtures',@catalogue; ...
 'Multiplicity matching is invariant under a permutation',@permutation; ...
 'A duplicated class cannot replace a missing class',@duplicate; ...
 'Different algebraic counts cannot be matched',@counts; ...
 'Class matching remains scale-relative at tiny and large common scales',@scales; ...
 'NaN, Inf and matrix-shaped representative outputs are rejected',@malformed; ...
 'Different quaternion directions in the same class are compared correctly',@orientations; ...
 'An invalid solver flag is not promoted to a verified reference match',@invalidInfo; ...
 'A zero eigenvector matrix cannot pass through a zero residual',@zeroVectors; ...
 'Every eigenvector column is checked, not just the Frobenius aggregate',@smallBadColumn; ...
 'Original quaternion D orientation is retained in the eigenpair equation',@vectorOrientation; ...
 'A nondiagonal D is not accepted as an eigenvalue diagonal',@nonDiagonal; ...
 'Repeated eigenvectors fail a required full basis despite small residuals',@rankDeficient; ...
 'A genuinely defective fixture does not require a nonexistent full basis',@defective; ...
 'The two-output eigenvalue multiset must agree with the one-output multiset',@twoOutputMismatch; ...
 'All analytic fixtures agree with the unchanged qstab right-spectrum solver',@primary; ...
 'Not requesting a legacy solver remains explicit and does not rank timings',@notRequested; ...
 'A correct external diagonal solver retains spectra and oriented eigenpairs',@correctLegacy; ...
 'Legacy spectral failures leave the successful qstab results intact',@wrongLegacy; ...
 'Legacy exceptions are retained without discarding other case data',@throwingLegacy; ...
 'Eigenvector-only errors retain a separately correct one-output spectrum',@vectorError; ...
 'Timing always requests exactly one output and checks every returned multiset',@timingContract; ...
 'Too-short eigen timing batches suppress speed ratios',@shortTiming; ...
 'A failed warm-up cannot disappear from a timing comparison',@failedTiming; ...
 'Case selection and invalid options fail before external solver execution',@selection; ...
 'Comparison preserves inputs, path, directory and the random state',@noSideEffects; ...
 'Comparison outputs contain no handles or captured private workspaces',@ordinary; ...
 'Full eigen-comparison data round-trip through the unchanged two-entry ZIP',@()roundtrip(writer)};
end
function F=fixture(id)
C=qstab_eigen_cases();F=C(strcmp({C.id},id));assert(isscalar(F));
end
function E=go(varargin)
E=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false,varargin{:});
end
function P=comp(a,b,c,d),P=struct('a',a,'b',b,'c',c,'d',d);end
function A=assess(F,E,varargin)
A=qstab_eigen_assess(F.matrix,F.reference,E,'Diagonalizable',F.diagonalizable,varargin{:});
end
function catalogue()
C=qstab_eigen_cases();assert(numel(C)==14&&numel(unique({C.id}))==14);
assert(isequaln(C,qstab_eigen_cases()));
for k=1:numel(C)
 n=size(C(k).matrix.a,1);assert(numel(C(k).reference)==n&&size(C(k).matrix.a,2)==n);
 assert(all(isfinite(C(k).reference))&&all(imag(C(k).reference)>=0));
end
end
function permutation()
y=[-2;-1+1i;-1+2i;-1+2i];x=y([3 1 4 2]);M=qstab_eigen_match(x,y);
assert(M.ok&&M.error==0&&numel(unique(M.permutation))==4);
assert(isequal(x(M.permutation),y));
end
function duplicate()
M=qstab_eigen_match([-1+1i;-1+1i;-1+3i],[-1+1i;-1+2i;-1+3i]);
assert(~M.ok&&strcmp(M.status,'MISMATCH')&&M.error>0);
end
function counts()
M=qstab_eigen_match([1;1],[1]);assert(~M.ok&&isinf(M.error));
M=qstab_eigen_match([],[]);assert(M.ok&&M.error==0);
end
function scales()
x=[-1+1i;-2+3i];y=x([2 1]);
for s=[1e-150,1,1e150]
 M=qstab_eigen_match(s*x,s*y);assert(M.ok&&M.error==0);
 N=qstab_eigen_match(s*[0;-2+3i],s*x);assert(~N.ok&&N.error>.1);
end
end
function malformed()
for x={NaN,Inf,ones(2),struct()}
 M=qstab_eigen_match(x{1},[1;2]);assert(~M.ok&&strcmp(M.status,'INVALID'));
end
fails(@()qstab_eigen_match(1,NaN));fails(@()qstab_eigen_match(1,1,0));
F=fixture('equal_real');R=assess(F,ones(3));assert(~R.ok&&strcmp(R.status,'INVALID'));
end
function orientations()
F=fixture('different_axes');D=F.matrix;
p=comp(diag(D.a),diag(D.b),diag(D.c),diag(D.d));R=assess(F,p);
assert(R.ok&&~R.standardForm&&R.spectrum.error==0);
end
function invalidInfo()
F=fixture('equal_real');R=assess(F,F.reference,'Info',struct('isValid',false));
assert(~R.ok&&strcmp(R.status,'UNRESOLVED')&&R.spectrum.ok);
end
function zeroVectors()
F=fixture('equal_real');R=assess(F,F.reference,'CheckVectors',true, ...
 'Vectors',zeros(3),'Diagonal',F.matrix);assert(~R.ok&&~R.vectors.allColumnsNonzero);
end
function smallBadColumn()
F=fixture('equal_real');V=eye(3);V(:,3)=1e-120*[1;0;0];
R=assess(F,F.reference,'CheckVectors',true,'Vectors',V,'Diagonal',F.matrix);
assert(~R.ok&&R.vectors.columnResiduals(3)>1e-3);
end
function vectorOrientation()
F=fixture('different_axes');D=F.matrix;P=comp(diag(D.a),diag(D.b),diag(D.c),diag(D.d));
R=assess(F,P,'CheckVectors',true,'Vectors',eye(3),'Diagonal',D);
assert(R.ok&&R.vectors.globalResidual==0&&R.vectors.fullBasisObserved);
end
function nonDiagonal()
F=fixture('equal_real');D=F.matrix;D.a(1,2)=1;
R=assess(F,F.reference,'CheckVectors',true,'Vectors',eye(3),'Diagonal',D);
assert(~R.ok&&strcmp(R.vectors.status,'INVALID'));
end
function rankDeficient()
F=fixture('repeated_real');V=repmat([1;0;0;0],1,4);
R=assess(F,F.reference,'CheckVectors',true,'Vectors',V,'Diagonal',-eye(4));
assert(~R.ok&&R.vectors.globalResidual==0&&~R.vectors.fullBasisObserved);
end
function defective()
F=fixture('jordan');V=repmat([1;0;0],1,3);
R=assess(F,F.reference,'CheckVectors',true,'Vectors',V,'Diagonal',-eye(3));
assert(R.ok&&~R.vectors.fullBasisObserved&&~R.vectors.fullBasisRequired);
end
function twoOutputMismatch()
F=fixture('equal_real');D=F.matrix;D.b(2,2)=1;
R=assess(F,F.reference,'CheckVectors',true,'Vectors',eye(3),'Diagonal',D);
assert(~R.ok&&R.spectrum.ok&&~R.vectors.diagonalAgreement.ok);
end
function primary()
E=go();assert(E.ok&&E.qstabOK&&E.matchedPrimary==14&&numel(E.cases)==14);
assert(~E.legacyRequested&&E.completed&&~E.allClassMultisetsAgree);
end
function notRequested()
E=go('Cases','equal_real','Timing',true);
assert(E.ok&&strcmp(E.cases.legacy.status,'NOT_REQUESTED'));
assert(strcmp(E.cases.timing.status,'NOT_APPLICABLE_NO_LEGACY')&&isnan(E.cases.timing.secondOverFirst));
end
function correctLegacy()
E=go('Cases',{'equal_real','different_axes','repeated_real'},'LegacyFunction',@diagonalLegacy);
assert(E.ok&&E.legacyOK&&E.allClassMultisetsAgree&&E.matchedLegacy==3);
for k=1:3,assert(E.cases(k).legacy.vectorAssessment.vectors.ok);end
end
function wrongLegacy()
E=go('Cases','equal_real','LegacyFunction',@wrongValues,'Timing',true,'CheckLegacyVectors',false);
assert(~E.ok&&E.qstabOK&&E.legacyIssues==1&&~E.allClassMultisetsAgree);
assert(strcmp(E.status,'COMPLETED_WITH_DIFFERENCES'));
assert(strcmp(E.cases.timing.status,'NOT_RUN_UNVERIFIED_SPECTRUM'));
end
function throwingLegacy()
E=go('Cases',{'equal_real','zero'},'LegacyFunction',@throws);
assert(~E.ok&&E.qstabOK&&E.matchedPrimary==2&&E.legacyIssues==2);
for k=1:2,assert(strcmp(E.cases(k).legacy.status,'ERROR'));end
end
function vectorError()
E=go('Cases','equal_real','LegacyFunction',@oneOutputOnly);
assert(~E.ok&&E.qstabOK&&E.cases.legacy.spectrumAssessment.ok);
assert(strcmp(E.cases.legacy.status,'EIGENVECTOR_ISSUE'));
end
function timingContract()
F=fixture('equal_real');T=qstab_eigen_timing(@onlyOne,@onlyOne,F.matrix,F.reference, ...
 'Repetitions',2,'WarmupRuns',1,'MinSampleTime',0,'MaxBatchCalls',1);
assert(T.ok&&numel(T.samples)==2&&T.firstCalls==4&&T.secondCalls==4);
assert(isequal(T.samples{1}.order,{'first','second'})&&isequal(T.samples{2}.order,{'second','first'}));
for k=1:2,assert(T.samples{k}.first.ok&&T.samples{k}.second.ok);end
end
function shortTiming()
F=fixture('equal_real');T=qstab_eigen_timing(@onlyOne,@onlyOne,F.matrix,F.reference, ...
 'Repetitions',1,'WarmupRuns',0,'MinSampleTime',1e6,'MaxBatchCalls',1);
assert(T.ok&&~T.eligible&&isnan(T.secondOverFirst)&&strcmp(T.status,'SHORT_BATCH'));
end
function failedTiming()
F=fixture('equal_real');T=qstab_eigen_timing(@onlyOne,@wrongValues,F.matrix,F.reference);
assert(~T.ok&&~T.eligible&&isnan(T.secondOverFirst)&&isempty(T.samples));
assert(numel(T.warmup)==1&&~T.warmup{1}.second.ok);
end
function selection()
E=go('Cases',{'scalar','equal_real'});assert(isequal({E.cases.id},{'equal_real','scalar'}));
fails(@()go('Cases','unknown','LegacyFunction',@throws));
fails(@()go('Cases',{'scalar','scalar'}));fails(@()go('Cases',{}));
fails(@()go('LegacyFunction','qeig'));fails(@()go('LegacyFunction',42));
fails(@()go('Repetitions',0));fails(@()go('Tolerance',NaN));
end
function noSideEffects()
C=qstab_eigen_cases();p=path;d=pwd;r=rng;
E=go('Cases','equal_real','LegacyFunction',@diagonalLegacy);
assert(E.ok&&strcmp(path,p)&&strcmp(pwd,d)&&isequaln(rng,r)&&isequaln(C,qstab_eigen_cases()));
assert(~any(structfun(@(x)x,E.environmentChangedByCalls)));
end
function ordinary()
E=go('Cases','equal_real','LegacyFunction',@diagonalLegacy);
assert(ischar(E.options.LegacyFunction)&&~isfield(E.legacySourceBefore,'workspace'));
assert(plain(E));
end
function roundtrip(writer)
base=tempname;mkdir(base);tidy=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
c=struct('group','fixture','name','Artificial eigen comparison persistence fixture', ...
 'status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
 'environment',struct('version','artificial','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','20000101T000000','runId','eigen-persistence-fixture', ...
 'nativeMATLABRun',false,'nativeQuaternionOK',false,'ok',true,'status','OVERALL_OK', ...
 'passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',c);
r.eigenComparison=go('Cases','equal_real','LegacyFunction',@diagonalLegacy);
r=writer(r,base,'qstab_test_results');folder=fullfile(base,'unpacked');mkdir(folder);unzip(r.archiveFile,folder);
f=dir(folder);f=f(~ismember({f.name},{'.','..'}));
assert(numel(f)==2&&isequal(sort({f.name}),sort({'report.mat','report.txt'})));
s=load(fullfile(folder,'report.mat'),'report');assert(isequaln(s.report.eigenComparison,r.eigenComparison));
end
function [V,D]=diagonalLegacy(A)
[a,b,c,d]=qstab.components(A);n=size(a,1);
if nargout<=1,V=comp(diag(a),diag(b),diag(c),diag(d));
else,V=eye(n);D=A;end
end
function E=onlyOne(A)
assert(nargout==1);E=diagonalLegacy(A);
end
function E=wrongValues(A)
E=diagonalLegacy(A);E.b(2)=E.b(1);
end
function [V,D]=oneOutputOnly(A)
if nargout>1,error('qstab:fixture:NoVectors','Synthetic two-output error.');end
V=diagonalLegacy(A);D=[];
end
function E=throws(~)
E=[];error('qstab:fixture:Deliberate','Synthetic external failure.'); %#ok<NASGU>
end
function ok=plain(x)
if isstruct(x)
 ok=true;f=fieldnames(x);for k=1:numel(x),for j=1:numel(f),ok=ok&&plain(x(k).(f{j}));end;end
elseif iscell(x),ok=all(cellfun(@plain,x(:)));
else,ok=isnumeric(x)||islogical(x)||ischar(x);end
end
function fails(f)
ok=false;try,f();catch,ok=true;end;assert(ok);
end
