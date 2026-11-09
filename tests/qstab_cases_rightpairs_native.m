function T=qstab_cases_rightpairs_native()
% Native MathWorks quaternion checks, explicitly separate from component tests.
T={ ...
 'Native right eigenpairs satisfy all fourteen analytic reference cases',@allCases; ...
 'Native and component eigenvectors retain the same right-eigenpair data',@sameData; ...
 'Explicit native output also accepts numeric matrix inputs',@numericToNative; ...
 'Native repeated real eigenvalues return a full quaternionic basis',@repeated; ...
 'Native canonicalization outside the i slice preserves oriented eigenpairs',@orientation; ...
 'Native primary-vector comparison remains distinct from an explicit external solver',@comparison};
end
function allCases()
C=qstab_eigen_cases();
for k=1:numel(C)
 F=C(k);P=F.matrix;A=quaternion(P.a,P.b,P.c,P.d);
 [V,D,I]=qstab.eigpairs(A);
 S=qstab_eigen_assess(A,F.reference,I.standardEigenvalues,'Info',I, ...
  'CheckVectors',true,'Vectors',V,'Diagonal',D,'Diagonalizable',F.diagonalizable);
 assert(S.ok,S.message);assert(isa(V,'quaternion')&&isa(D,'quaternion'));
end
end
function sameData()
A=quaternion([-1 .2;0 -2],[1 0;0 0],[0 .1;0 2],[0 0;0 0]);
[V,D,I]=qstab.eigpairs(A);[W,E,J]=qstab.eigpairs(A,'Output','components');
[a,b,c,d]=parts(V);assert(isequaln(W,struct('a',a,'b',b,'c',c,'d',d)));
[a,b,c,d]=parts(D);assert(isequaln(E,struct('a',a,'b',b,'c',c,'d',d)));
I.outputType='components';assert(isequaln(I,J));
end
function numericToNative()
[V,D,I]=qstab.eigpairs(-eye(3),'Output','quaternion');
assert(isa(V,'quaternion')&&isa(D,'quaternion')&&I.fullBasisObserved);
assert(I.quaternionBasisRank==3);
end
function repeated()
A=quaternion(-eye(4),zeros(4),zeros(4),zeros(4));[V,D,I]=qstab.eigpairs(A);
assert(I.fullBasisObserved&&I.quaternionBasisRank==4);
assert(norm(qstab.chi(A)*qstab.chi(V)-qstab.chi(V)*qstab.chi(D),'fro')==0);
end
function orientation()
A=quaternion(-2,1,-2,.5);[V,D,I]=qstab.eigpairs(A);
[a,b,c,d]=parts(D);assert(abs(a+2)<1e-13&&abs(b-sqrt(5.25))<1e-13&&c==0&&d==0);
assert(I.eigenpairsValid&&max(I.columnResiduals)<1e-13&&isa(V,'quaternion'));
end
function comparison()
E=run_qstab_eigen_comparison('Cases',{'equal_real','repeated_real'}, ...
 'CheckPrimaryVectors',true,'UseQuaternion',true,'Verbose',false,'LegacyFunction',@synthetic);
assert(E.ok&&E.qstabOK&&E.primaryEigenpairsOK&&E.legacyOK&&E.allClassMultisetsAgree);
end
function [V,D]=synthetic(A)
% Synthetic correct interface, not the user's private qeig.
[V,D,I]=qstab.eigpairs(A);if nargout<=1
 e=I.standardEigenvalues;V=quaternion(real(e),imag(e),zeros(size(e)),zeros(size(e)));
end
end
