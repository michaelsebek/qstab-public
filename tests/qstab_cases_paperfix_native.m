function T=qstab_cases_paperfix_native()
%QSTAB_CASES_PAPERFIX_NATIVE Source evidence with genuine MathWorks quaternion.
T={ ...
 'Native Schur source witness retains separate unresolved and spectral evidence',@schurWitness; ...
 'Native correlated source certificates preserve the same shared input family',@correlated};
end
function E=go(id,native)
E=run_qstab_papers('Cases',id,'UseQuaternion',native,'Verbose',false,'MaxTime',Inf);
assert(E.ok&&E.fullyMatched);
end
function schurWitness()
E=go('qb03_schur_tetrahedron',true);D=E.cases.output;
assert(E.nativeQuaternionInputs&&strcmp(D.result.status,'undetermined'));
assert(strcmp(D.centreResult.status,'notStable')&&D.centreResult.spectrumComputed);
assert(isequal(D.centre,D.centreFromWeights)&&sum(D.centreWeights)==1);
B=go('qb03_schur_tetrahedron',false);
a=sort(real(D.centreResult.diagnostics.representatives(:)));
b=sort(real(B.cases.output.centreResult.diagnostics.representatives(:)));
assert(norm(a-b)<1e-9&&max(a)>1);
end
function correlated()
A=go('qb03_correlation',true);B=go('qb03_correlation',false);
assert(A.nativeQuaternionInputs);a=A.cases.output;b=B.cases.output;
assert(strcmp(a.correlated.status,b.correlated.status));
L=a.correlatedLyapunov;M=b.correlatedLyapunov;
assert(strcmp(L.status,'stable')&&L.certificateComplete&&L.unresolvedCells==0);
assert(L.certifiedCells==M.certifiedCells&&numel(L.certificates)==L.certifiedCells);
assert(abs(L.certificate.minimumSlack-M.certificate.minimumSlack)<1e-8);
assert(isequal(a.independentNorm,b.independentNorm)&&strcmp(a.independent.status,'notStable'));
end
