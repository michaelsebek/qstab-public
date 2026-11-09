function T=qstab_cases_paperfix(writer)
%QSTAB_CASES_PAPERFIX Explicit evidence for source cases, not changed solvers.
T={ ...
 'Schur equality remains unresolved in the original root-free nominal test',@schurEquality; ...
 'The Schur source centre is an admissible spectral counterexample',@schurCentre; ...
 'The original Schur zero-exclusion result is retained without relabelling',@schurRetained; ...
 'The correlated source segment has stored whole-family Lyapunov coverage',@correlatedCoverage; ...
 'The original correlated zero-exclusion output is retained independently',@correlatedRetained; ...
 'The independent unstable member is outside the shared-parameter segment',@independentMember; ...
 'An insufficient positive-evidence budget remains explicitly unresolved',@limitedBudget; ...
 'An unresolved paper check identifies its actual and reference statuses',@failureDetail; ...
 'Alternative evidence does not modify source coefficients or references',@sourcePreserved; ...
 'Both repaired source cases round-trip with their separate method outputs',@()roundtrip(writer)};
end
function C=source(id)
C=qstab_paper_cases();C=C(strcmp({C.id},id));assert(isscalar(C));
end
function E=go(id,varargin)
E=run_qstab_papers('Cases',id,'UseQuaternion',false,'Verbose',false,'MaxTime',Inf,varargin{:});
end
function P=poly(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function V=polys(C)
V=cell(size(C));for j=1:numel(C),V{j}=poly(C{j});end
end
function schurEquality()
C=source('qb03_schur_tetrahedron');M=mean(cat(3,C.data.vertices{:}),3);
N=qstab.normpoly(poly(M));assert(isequal(N,[1 8 14 -8 1]));
[~,R] = qstab.isstable(poly(M),'Variable','z');
assert(strcmp(R.status,'undetermined')&&strcmp(R.method,'schur')&&~R.spectrumComputed);
assert(R.diagnostics.test.reflectionGaps(1)==0);
end
function schurCentre()
E=go('qb03_schur_tetrahedron');assert(E.ok&&E.fullyMatched);D=E.cases.output;
assert(all(D.centreWeights>=0)&&sum(D.centreWeights)==1);
assert(isequal(D.centre,D.centreFromWeights)&&isequal(D.centre,E.catalogue.reference.centre));
S=D.centreResult;assert(strcmp(S.status,'notStable')&&S.spectrumComputed&&strcmp(S.method,'spectral'));
lambda=S.diagnostics.representatives(:);
assert(numel(lambda)==2&&max(abs(imag(lambda)))<1e-9);
assert(norm(sort(real(lambda))-[2-sqrt(5);2+sqrt(5)])<1e-8);
end
function schurRetained()
E=go('qb03_schur_tetrahedron');D=E.cases.output;
[~,R] = qstab.zeroexclude(polys(E.catalogue.data.vertices),'Variable','z','MaxTime',Inf,'MaxCells',8192);
assert(strcmp(D.result.status,R.status)&&strcmp(D.result.status,'undetermined'));
assert(strcmp(D.result.nominal.method,'schur')&&~D.result.nominal.spectrumComputed);
assert(D.result.processedCells==R.processedCells&&isempty(D.result.witness));
assert(strcmp(D.centreResult.status,'notStable')&&E.fullyMatched);
end
function correlatedCoverage()
E=go('qb03_correlation');assert(E.ok&&E.fullyMatched);D=E.cases.output;R=D.correlatedLyapunov;
assert(strcmp(R.status,'stable')&&strcmp(R.method,'subdivisionLyapunovPerturbation'));
assert(R.certificateComplete&&R.unresolvedCells==0&&R.remainingCells==0);
assert(R.certifiedCells>0&&numel(R.certificates)==R.certifiedCells);
assert(R.options.StoreCertificates&&R.options.MaxCells==8192&&isinf(R.options.MaxTime));
V=E.catalogue.data.vertices;intervals=zeros(numel(R.certificates),2);
for j=1:numel(R.certificates)
 A=R.certificates{j};W=A.cell.weights;
 assert(size(W,1)==2&&all(W(:)>=0)&&all(sum(W,1)==1));
 P=zeros(size(A.cell.vertices));
 for col=1:size(W,2)
  for k=1:2,P(:,:,col)=P(:,:,col)+W(k,col)*V{k};end
 end
 assert(isequal(P,A.cell.vertices)&&A.accepted&&A.slack>0);
 theta=[-1 1]*W;intervals(j,:)=[min(theta),max(theta)];
end
intervals=sortrows(intervals,1);
assert(intervals(1,1)==-1&&intervals(end,2)==1&&all(intervals(:,2)>intervals(:,1)));
assert(all(intervals(1:end-1,2)==intervals(2:end,1)));
end
function correlatedRetained()
E=go('qb03_correlation');D=E.cases.output;
[~,R] = qstab.zeroexclude(polys(E.catalogue.data.vertices),'MaxTime',Inf,'MaxCells',8192);
assert(strcmp(D.correlated.status,R.status)&&D.correlated.processedCells==R.processedCells);
assert(strcmp(D.correlated.method,'normValueZeroExclusion'));
assert(isequaln(D.correlated.unresolvedRegions,R.unresolvedRegions));
assert(any(strcmp(D.correlated.status,{'stable','undetermined'})));
assert(strcmp(D.correlatedLyapunov.status,'stable'));
end
function independentMember()
E=go('qb03_correlation');D=E.cases.output;C=E.catalogue.data;
V=cat(3,C.vertices{:});L=min(V,[],3);U=max(V,[],3);P=C.independent;
assert(all(P(:)>=L(:))&&all(P(:)<=U(:)));
assert(P(3,2)~=3*P(3,1)); % The two occurrences cannot share one real parameter.
assert(strcmp(D.independent.status,'notStable')&&isequal(D.independentNorm,[3 -2 17 2 1]));
end
function limitedBudget()
E=go('qb03_correlation','MaxCells',1);D=E.cases.output;
assert(E.ok&&~E.fullyMatched&&E.unresolved==1&&strcmp(E.cases.status,'UNRESOLVED'));
assert(strcmp(D.correlatedLyapunov.status,'undetermined')&&~D.correlatedLyapunov.certificateComplete);
assert(D.correlatedLyapunov.options.MaxCells==1&&D.correlated.options.MaxCells==1);
assert(D.correlatedLyapunov.unresolvedCells>0&&~isempty(D.correlatedLyapunov.unresolvedRegions));
assert(strcmp(E.cases.checks{1}.status,'UNRESOLVED'));
end
function failureDetail()
E=go('qb03_correlation','MaxCells',1);m=E.cases.message;
assert(~isempty(strfind(m,'UNRESOLVED'))&&~isempty(strfind(m,'actual=undetermined')));
assert(~isempty(strfind(m,'reference=stable'))&&~isempty(strfind(m,'Check 1')));
end
function sourcePreserved()
C=qstab_paper_cases();ids={'qb03_schur_tetrahedron','qb03_correlation'};
E=go(ids,'MaxCells',1);after=qstab_paper_cases();assert(isequaln(C,after));
for j=1:numel(E.catalogue)
 k=find(strcmp({C.id},E.catalogue(j).id));assert(isequaln(C(k),E.catalogue(j)));
end
assert(~E.fullyMatched&&E.ok);
end
function roundtrip(writer)
base=tempname;mkdir(base);cleanup=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
item=struct('group','fixture','name','Artificial source-evidence persistence fixture', ...
 'status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
 'environment',struct('version','artificial','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','20000101T000000','runId','paperfix-persistence-fixture', ...
 'nativeMATLABRun',false,'nativeQuaternionOK',false,'ok',true,'status','OVERALL_OK', ...
 'passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',item);
r.paperExamples=go({'qb03_schur_tetrahedron','qb03_correlation'});
assert(r.paperExamples.ok&&r.paperExamples.fullyMatched);
r=writer(r,base,'qstab_test_results');folder=fullfile(base,'unpacked');mkdir(folder);unzip(r.archiveFile,folder);
f=dir(folder);f=f(~ismember({f.name},{'.','..'}));
assert(numel(f)==2&&isequal(sort({f.name}),sort({'report.mat','report.txt'})));
s=load(fullfile(folder,'report.mat'),'report');assert(isequaln(s.report.paperExamples,r.paperExamples));
assert(strcmp(r.archiveFile,fullfile(base,'qstab_test_results.zip')));
end
