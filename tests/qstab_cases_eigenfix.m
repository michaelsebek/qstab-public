function T=qstab_cases_eigenfix(writer)
%QSTAB_CASES_EIGENFIX Preserve a fixed outer schema through every case path.
% The synthetic external functions below are test fixtures, NOT old qeig.
T={ ...
 'The first comparison row is assigned into a matching six-field schema',@firstRow; ...
 'A selected final case fills exactly one row with no preallocation placeholders',@selectedFinal; ...
 'An external error in the first row cannot discard a later successful row',@errorThenSuccess; ...
 'A successful row followed by an error retains both distinct nested payloads',@successThenError; ...
 'An eigenvector-only failure leaves the next row and one-output spectrum intact',@vectorFailure; ...
 'Timed and untimed error branches fit the same preallocated case schema',@timingBranches; ...
 'Invalid external output remains a diagnostic row rather than an assignment error',@invalidOutput; ...
 'Mixed successful and failed comparison rows survive the fixed two-entry ZIP',@()roundtrip(writer)};
end
function E=go(varargin)
E=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false,varargin{:});
end
function assertSchema(E,ids)
names={'id';'primary';'legacy';'mutualMatch';'classMultisetsAgree';'timing'};
assert(isequal(size(E.cases),[1 numel(ids)]));
assert(isequal(fieldnames(E.cases),names));
assert(isequal({E.cases.id},ids));
assert(E.completed);
for k=1:numel(E.cases)
 assert(isequal(fieldnames(E.cases(k)),names));
 assert(isstruct(E.cases(k).primary)&&isscalar(E.cases(k).primary));
 assert(isstruct(E.cases(k).legacy)&&isscalar(E.cases(k).legacy));
 assert(isstruct(E.cases(k).timing)&&isscalar(E.cases(k).timing));
end
end
function firstRow()
E=go('Cases','equal_real');assertSchema(E,{'equal_real'});
assert(E.ok&&E.qstabOK&&E.matchedPrimary==1);
assert(strcmp(E.cases.legacy.status,'NOT_REQUESTED'));
assert(strcmp(E.cases.timing.status,'NOT_REQUESTED'));
end
function selectedFinal()
E=go('Cases','dense16');assertSchema(E,{'dense16'});
assert(E.ok&&E.matchedPrimary==1&&numel(E.catalogue)==1);
assert(numel(E.cases.primary.canonical)==16);
end
function errorThenSuccess()
E=go('Cases',{'zero','equal_real'},'LegacyFunction',@errorOnNonzero);
assertSchema(E,{'equal_real','zero'});
assert(~E.ok&&E.qstabOK&&E.matchedPrimary==2&&E.matchedLegacy==1&&E.legacyIssues==1);
assert(strcmp(E.status,'COMPLETED_WITH_DIFFERENCES'));
assert(strcmp(E.cases(1).legacy.status,'ERROR'));
assert(strcmp(E.cases(1).legacy.identifier,'qstab:fixture:EigenFixNonzero'));
assert(isempty(E.cases(1).legacy.spectrumAssessment));
assert(E.cases(2).legacy.ok&&E.cases(2).classMultisetsAgree);
assert(~isfield(E.cases(2).legacy,'identifier'));
assert(isempty(E.cases(1).mutualMatch)&&E.cases(2).mutualMatch.ok);
end
function successThenError()
E=go('Cases',{'zero','equal_real'},'LegacyFunction',@errorOnZero);
assertSchema(E,{'equal_real','zero'});
assert(~E.ok&&E.qstabOK&&E.matchedPrimary==2&&E.matchedLegacy==1&&E.legacyIssues==1);
assert(E.cases(1).legacy.ok&&E.cases(1).classMultisetsAgree);
assert(strcmp(E.cases(2).legacy.status,'ERROR'));
assert(strcmp(E.cases(2).legacy.identifier,'qstab:fixture:EigenFixZero'));
assert(isempty(E.cases(2).legacy.spectrumAssessment)&&isempty(E.cases(2).mutualMatch));
end
function vectorFailure()
E=go('Cases',{'equal_real','zero'},'LegacyFunction',@vectorErrorOnNonzero);
assertSchema(E,{'equal_real','zero'});
assert(~E.ok&&E.qstabOK&&E.allClassMultisetsAgree&&E.matchedLegacy==1);
assert(E.cases(1).legacy.spectrumAssessment.ok);
assert(strcmp(E.cases(1).legacy.status,'EIGENVECTOR_ISSUE'));
assert(strcmp(E.cases(1).legacy.vectorAssessment.identifier,'qstab:fixture:EigenFixVectors'));
assert(E.cases(2).legacy.ok&&E.cases(2).legacy.vectorAssessment.vectors.ok);
end
function timingBranches()
E=go('Cases',{'equal_real','zero'},'LegacyFunction',@errorOnZero, ...
 'Timing',true,'Repetitions',1,'WarmupRuns',0,'MinSampleTime',1e6,'MaxBatchCalls',1);
assertSchema(E,{'equal_real','zero'});
assert(E.qstabOK&&~E.ok&&E.timingEligible==0);
assert(strcmp(E.cases(1).timing.status,'SHORT_BATCH'));
assert(E.cases(1).timing.ok&&~isempty(E.cases(1).timing.samples));
assert(strcmp(E.cases(2).timing.status,'NOT_RUN_UNVERIFIED_SPECTRUM'));
assert(isnan(E.cases(1).timing.secondOverFirst)&&isnan(E.cases(2).timing.secondOverFirst));
assert(~isfield(E.cases(2).timing,'samples'));
end
function invalidOutput()
E=go('Cases',{'equal_real','zero'},'LegacyFunction',@invalidOnNonzero, ...
 'CheckLegacyVectors',false);
assertSchema(E,{'equal_real','zero'});
assert(E.qstabOK&&~E.ok&&E.matchedLegacy==1&&E.legacyIssues==1);
assert(~E.cases(1).legacy.spectrumAssessment.ok);
assert(strcmp(E.cases(1).legacy.status,'INVALID'));
assert(E.cases(2).legacy.ok&&E.cases(2).classMultisetsAgree);
end
function roundtrip(writer)
base=tempname;mkdir(base);tidy=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
c=struct('group','fixture','name','Artificial mixed-schema persistence fixture', ...
 'status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
 'environment',struct('version','artificial','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','20000101T000000','runId','eigen-schema-persistence-fixture', ...
 'nativeMATLABRun',false,'nativeQuaternionOK',false,'ok',false,'status','EIGEN_COMPARISON_ISSUES', ...
 'passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',c);
r.eigenComparison=go('Cases',{'equal_real','zero'},'LegacyFunction',@errorOnNonzero);
r=writer(r,base,'qstab_test_results');folder=fullfile(base,'unpacked');mkdir(folder);
unzip(r.archiveFile,folder);f=dir(folder);f=f(~ismember({f.name},{'.','..'}));
assert(numel(f)==2&&isequal(sort({f.name}),sort({'report.mat','report.txt'})));
s=load(fullfile(folder,'report.mat'),'report');
assert(isequaln(s.report.eigenComparison,r.eigenComparison));
assert(strcmp(s.report.status,'EIGEN_COMPARISON_ISSUES'));
assertSchema(s.report.eigenComparison,{'equal_real','zero'});
end
function [V,D]=diagonalFixture(A)
[a,b,c,d]=qstab.components(A);n=size(a,1);
if nargout<=1,V=struct('a',diag(a),'b',diag(b),'c',diag(c),'d',diag(d));D=[];
else,V=eye(n);D=A;end
end
function tf=isZero(A)
[a,b,c,d]=qstab.components(A);tf=~any([a(:);b(:);c(:);d(:)]);
end
function [V,D]=errorOnNonzero(A)
if ~isZero(A),error('qstab:fixture:EigenFixNonzero','Synthetic first-case error.');end
if nargout<=1,V=diagonalFixture(A);D=[];else,[V,D]=diagonalFixture(A);end
end
function [V,D]=errorOnZero(A)
if isZero(A),error('qstab:fixture:EigenFixZero','Synthetic later-case error.');end
if nargout<=1,V=diagonalFixture(A);D=[];else,[V,D]=diagonalFixture(A);end
end
function [V,D]=vectorErrorOnNonzero(A)
if nargout>1&&~isZero(A)
 error('qstab:fixture:EigenFixVectors','Synthetic two-output error.');
end
if nargout<=1,V=diagonalFixture(A);D=[];else,[V,D]=diagonalFixture(A);end
end
function V=invalidOnNonzero(A)
if isZero(A),V=diagonalFixture(A);else,V=NaN(size(A.a,1),1);end
end
