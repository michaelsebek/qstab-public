function T=qstab_cases_eigencomparison_native()
T={ ...
 'Native eigen reference fixtures agree with their explicit component forms',@fixtures; ...
 'Native oriented diagonal eigenvectors are checked without canonicalizing D',@oriented; ...
 'Native equal-real-part repeated classes retain algebraic multiplicity',@multiplicity; ...
 'Native comparison results remain ordinary serializable data',@serialization};
end
function E=go(varargin)
E=run_qstab_eigen_comparison('UseQuaternion',true,'Verbose',false,varargin{:});
end
function fixtures()
A=go();B=run_qstab_eigen_comparison('UseQuaternion',false,'Verbose',false);
assert(A.ok&&B.ok&&A.nativeQuaternionInputs&&A.matchedPrimary==14);
for k=1:14
 M=qstab_eigen_match(A.cases(k).primary.canonical,B.cases(k).primary.canonical);
 assert(M.ok);
end
end
function oriented()
E=go('Cases','different_axes','LegacyFunction',@nativeDiagonal);
assert(E.ok&&E.allClassMultisetsAgree&&E.cases.legacy.vectorAssessment.vectors.globalResidual==0);
assert(~E.cases.legacy.spectrumAssessment.standardForm);
end
function multiplicity()
E=go('Cases',{'equal_real','repeated_real','repeated_nonreal','mixed_multiplicity'}, ...
 'LegacyFunction',@nativeDiagonal);
assert(E.ok&&E.matchedLegacy==4&&E.allClassMultisetsAgree);
end
function serialization()
E=go('Cases','equal_real','LegacyFunction',@nativeDiagonal);
assert(E.ok&&isstruct(E.cases.legacy.spectrumAssessment.components));
f=[tempname '.mat'];tidy=onCleanup(@()delete(f)); %#ok<NASGU>
save(f,'E','-v7');R=load(f,'E');assert(isequaln(R.E,E));
end
function [V,D]=nativeDiagonal(A)
[a,b,c,d]=qstab.components(A);n=size(a,1);
if nargout<=1,V=quaternion(diag(a),diag(b),diag(c),diag(d));
else,V=quaternion(eye(n),zeros(n),zeros(n),zeros(n));D=A;end
end
