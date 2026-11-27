function cases=qstab_cases_oneparameter_examples()
%QSTAB_CASES_ONEPARAMETER_EXAMPLES Article-data and reproduction integration.
cases={...
 'Scalar article catalogue preserves all nine frequency and fourteen ray fixtures',@catalogue;...
 'Rank-one article example runs through the public result-first functions',@rankOne;...
 'Nearest-gain sphere example checks original polynomial witnesses',@sphereCase;...
 'Tangencies remain explicitly unresolved unless a proof is available',@tangencies;...
 'One-parameter runner preserves caller state and rejects unknown selections',@stateAndOptions};
end
function catalogue()
C=qstab_oneparameter_cases();assert(numel(C)==13);
assert(isequal({C(1:7).id},{'E1','E2','E3','E4','E5','E6','E7'}));
np=0;nr=0;for k=1:numel(C),if isfield(C(k).data,'problems'),np=np+numel(C(k).data.problems);end;nr=nr+numel(C(k).rays);end
assert(np==9&&nr==14);
assert(isequal(C(1).rays(1).direction(:),[0;1;0;1]));
assert(strcmp(C(7).data.sourceId,'K')&&C(7).data.frequencyScale==1e6);
assert(C(strcmp({C.id},'E5V')).rays.reference==1);
end
function rankOne()
R=runCase('E2');assert(R.ok&&~R.articleComplete&&~R.exact.executed);
D=R.cases.output;assert(abs(D.problems{1}.frequencyRadius-sqrt(5/8))<1e-10);
assert(abs(D.rays{1}.interval(2)-.25)<1e-8);
assert(norm(D.originalResidual)<1e-10&&sum(D.perturbationSingularValues>1e-10)==1);
assert(strcmp(R.exact.status,'NOT_REQUESTED'));
end
function sphereCase()
R=runCase('E3');assert(R.ok);D=R.cases.output;
assert(size(D.sphereParameters,2)==7&&max(D.sphereResiduals)<1e-10);
assert(max(abs(sqrt(sum(D.sphereParameters.^2,1))-sqrt(3)/2))<1e-12);
assert(~R.articleComplete&&~R.exact.executed);
end
function tangencies()
R=runCase({'E5','E5V','E5C'});assert(R.ok&&~R.articleComplete);
for k=1:numel(R.cases)
 D=R.cases(k).output;I=D.rays{1}.info;
 if strcmp(I.positiveStatus,'computed'),assert(abs(D.rays{1}.interval(2)-1)<2e-7);
 else,assert(isnan(D.rays{1}.interval(2))&&D.unresolvedRays>=1);end
 assert(~D.rays{1}.exact.executed);
end
end
function stateAndOptions()
p=path;folder=pwd;r=rng;R=runCase('A01');assert(R.ok);
assert(strcmp(path,p)&&strcmp(pwd,folder)&&isequal(rng,r));
assert(~R.articleComplete&&numel(R.sourceHashes)==2);
bad=false;try,runCase('E99');catch e,bad=strcmp(e.identifier,'qstab:oneparameter:Selection');end
assert(bad,'An unknown article key must fail explicitly.');
end
function R=runCase(key)
R=run_qstab_oneparameter('Cases',key,'UseQuaternion',false,'ExactCertificates',false, ...
 'Save',false,'Plots',false,'Verbose',false,'MaxTime',.25,'MaxCells',64);
end
