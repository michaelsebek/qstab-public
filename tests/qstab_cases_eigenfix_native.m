function T=qstab_cases_eigenfix_native()
%QSTAB_CASES_EIGENFIX_NATIVE Real MathWorks quaternion inputs, synthetic legacy.
T={ ...
 'Native first-row comparison preserves the fixed case schema without a legacy solver',@firstRow; ...
 'Native successful and throwing legacy paths coexist in one comparison array',@mixedRows};
end
function firstRow()
E=run_qstab_eigen_comparison('UseQuaternion',true,'Verbose',false,'Cases','equal_real');
assert(E.completed&&E.ok&&E.qstabOK&&E.matchedPrimary==1);
assert(isequal(size(E.cases),[1 1]));
assert(isequal(fieldnames(E.cases), ...
 {'id';'primary';'legacy';'mutualMatch';'classMultisetsAgree';'timing'}));
assert(strcmp(E.cases.id,'equal_real'));
assert(strcmp(E.cases.legacy.status,'NOT_REQUESTED'));
end
function mixedRows()
E=run_qstab_eigen_comparison('UseQuaternion',true,'Verbose',false, ...
 'Cases',{'equal_real','zero'},'LegacyFunction',@nativeErrorOnZero);
assert(E.completed&&~E.ok&&E.qstabOK&&E.matchedPrimary==2);
assert(E.matchedLegacy==1&&E.legacyIssues==1);
assert(isequal({E.cases.id},{'equal_real','zero'}));
assert(E.cases(1).legacy.ok&&E.cases(1).legacy.vectorAssessment.vectors.ok);
assert(strcmp(E.cases(2).legacy.status,'ERROR'));
assert(strcmp(E.cases(2).legacy.identifier,'qstab:fixture:NativeEigenFixZero'));
assert(isempty(E.cases(2).mutualMatch));
end
function [V,D]=nativeErrorOnZero(A)
assert(strcmp(class(A),'quaternion'));
[a,b,c,d]=qstab.components(A);n=size(a,1);
if ~any([a(:);b(:);c(:);d(:)])
 error('qstab:fixture:NativeEigenFixZero','Synthetic native-input zero-matrix error.');
end
if nargout<=1,V=quaternion(diag(a),diag(b),diag(c),diag(d));D=[];
else,V=quaternion(eye(n),zeros(n),zeros(n),zeros(n));D=A;end
end
