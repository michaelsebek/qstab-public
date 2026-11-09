function cases = qstab_cases_routhsingular()
%QSTAB_CASES_ROUTHSINGULAR Root-free quartic decisions at singular pivots.
% The original unstable quadratic fixture is deliberately retained. These
% tests distinguish a singular Routh pivot from a stability-boundary root.
cases={ ...
 'A zero Routh pivot with negative third Hurwitz determinant is unstable',@centerWitness; ...
 'Auto, coefficient and Routh methods preserve the root-free witness',@methodAliases; ...
 'Common small, large and negative scales preserve the quartic witness',@commonScales; ...
 'Nonreal nonmonic common factors retain the same unstable classes',@commonQuaternionFactors; ...
 'Tiny pivot perturbations do not mask a well separated negative determinant',@tinyPivots; ...
 'Unresolved quartics and exact structural axis witnesses are distinguished',@boundaryBand; ...
 'Stable quadratic corners retain ordinary positive Routh decisions',@stableCorners; ...
 'Schur and nonquartic tests do not use the quartic Hurwitz fallback',@otherDomains; ...
 'Interval routes retain the same admissible unstable midpoint witness',@intervalRoutes; ...
 'Nine analytic zero-pivot quadratics agree with independent spectral checks',@analyticFamily};
end

function centerWitness()
C=center();[~,R] = qstab.isstable(poly(C));
unstable(R);D=R.diagnostics.test;W=D.hurwitzFallback;
assert(D.criticalRow==3 && D.firstColumn(end)==0);
assert(strcmp(W.name,'quarticThirdHurwitzDeterminant') && W.order==3);
assert(W.normalizedDeterminant==-.25 && W.negativeGap>0);
assert(W.provesRightHalfPlaneRoot && ~W.validated);
N=qstab.normpoly(poly(C));assert(isequal(N,[4,4,2,2,1]));
% Direct determinant, independent of the three-term implementation.
a=W.coefficientsDescending;
H=[a(2),a(4),0;a(1),a(3),a(5);0,a(2),a(4)];
assert(abs(det(H)-W.normalizedDeterminant)<1e-14);
end

function methodAliases()
C=center();P=poly(C);
for method={'auto','coefficient','routh'}
 [~,R] = qstab.isstable(P,'Method',method{1});unstable(R);
 assert(R.diagnostics.test.hurwitzFallback.provesRightHalfPlaneRoot);
end
end

function commonScales()
C=center();[~,baseline] = qstab.isstable(poly(C));
for s=[1e-140,-1e-140,1,-1,1e140,-1e140]
 [~,R] = qstab.isstable(poly(s*C));unstable(R);
 assert(abs(R.diagnostics.test.hurwitzFallback.normalizedDeterminant- ...
  baseline.diagnostics.test.hurwitzFallback.normalizedDeterminant)<1e-13);
end
end

function commonQuaternionFactors()
C=center();q=[1;2;3;4];L=zeros(size(C));R=L;
for k=1:size(C,2)
 L(:,k)=multiply(q,C(:,k));R(:,k)=multiply(C(:,k),q);
end
for P={poly(L),poly(R)}
 [~,S] = qstab.isstable(P{1});unstable(S);
 W=S.diagnostics.test.hurwitzFallback;
 assert(W.provesRightHalfPlaneRoot && W.normalizedDeterminant < -.24);
end
end

function tinyPivots()
for delta=[-1e-12,0,1e-12]
 C=center();C(2,2)=1+delta;[~,R] = qstab.isstable(poly(C));unstable(R);
 W=R.diagnostics.test.hurwitzFallback;
 assert(W.provesRightHalfPlaneRoot && W.negativeGap>.2);
end
end

function boundaryBand()
for delta=[-1e-12,0,1e-12]
 C=center();C(3,2)=1;C(4,2)=1+delta;
 [~,R] = qstab.isstable(poly(C),'Tolerance',1e-10);
 assert(strcmp(R.status,'undetermined') && isempty(R.isStable));
 assert(~R.spectrumComputed);
 W=R.diagnostics.test.hurwitzFallback;
 assert(~W.provesRightHalfPlaneRoot && W.negativeGap<=0);
end
[tf,R] = qstab.isstable([1,0,1]);
assert(tf==0 && strcmp(R.status,'notStable') && ~R.spectrumComputed && R.diagnostics.boundaryWitness.exactForSuppliedComponents);
assert(~R.diagnostics.test.hurwitzFallback.provesRightHalfPlaneRoot);
end

function stableCorners()
for x=[-2,2]
 for y=[-2,2]
  C=center();C(3,2)=x;C(4,2)=y;[~,R] = qstab.isstable(poly(C));
  assert(strcmp(R.status,'stable') && R.isStable && ~R.spectrumComputed);
  assert(all(R.diagnostics.test.firstColumn>0));
  assert(~isfield(R.diagnostics.test,'hurwitzFallback'));
 end
end
end

function otherDomains()
[tf,R] = qstab.isstable([1,0,1],'Variable','z');
assert(tf==0 && strcmp(R.status,'notStable') && ~R.spectrumComputed && R.diagnostics.boundaryWitness.exactForSuppliedComponents);
assert(~isfield(R.diagnostics.test,'hurwitzFallback'));
[~,R] = qstab.isstable([1,0,1,0,1],'Variable','s');
assert(strcmp(R.status,'undetermined') && ~R.spectrumComputed);
assert(~isfield(R.diagnostics.test,'hurwitzFallback'));
end

function intervalRoutes()
C=center();L=C;U=C;L(3:4,2)=-2;U(3:4,2)=2;V={poly(L),poly(U)};
for mode={'none','quadratic-faces','kharitonov-cube','kharitonov-faces'}
 [~,R] = qstab.zeroexclude(V,'Uncertainty','interval','Reduction',mode{1},'MaxTime',Inf);
 assert(strcmp(R.status,'notStable') && ~R.isStable);
 if strcmp(mode{1},'none'),actual=R.witness.components;
 else,actual=R.witness.originalComponents;end
 assert(isequal(actual,C));
 assert(~R.witness.nominal.spectrumComputed);
 assert(R.witness.nominal.diagnostics.test.hurwitzFallback.provesRightHalfPlaneRoot);
end
[~,R] = qstab.robuststable(V,'Uncertainty','interval','Method','zeroexclusion', ...
 'Reduction','quadratic-faces','MaxTime',Inf);
assert(strcmp(R.status,'notStable') && isequal(R.witness.originalComponents,C));
end

function analyticFamily()
% a=alpha+gamma*i; b=beta*i; beta=alpha*(alpha^2+gamma^2)/gamma.
% Thus Delta2=0 and Delta3=-4*alpha^2*beta^2<0 exactly.
for alpha=[.5,1,2]
 for gamma=[.5,1,2]
  beta=alpha*(alpha^2+gamma^2)/gamma;
  C=[0,alpha,1;beta,gamma,0;0,0,0;0,0,0];P=poly(C);
  [~,R] = qstab.isstable(P);unstable(R);
  [~,S] = qstab.isstable(P,'Method','spectral');
  assert(strcmp(S.status,'notStable') && S.spectrumComputed);
  assert(R.diagnostics.test.hurwitzFallback.provesRightHalfPlaneRoot);
 end
end
end

function C=center()
C=[0,1,1;2,1,0;0,0,0;0,0,0];
end
function unstable(R)
assert(strcmp(R.status,'notStable') && isequal(R.isStable,false),R.message);
assert(~R.spectrumComputed && strcmp(R.method,'routh'));
end
function P=poly(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
 'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function z=multiply(x,y)
% Hamilton product, written explicitly to keep this fixture independent.
z=[x(1)*y(1)-x(2)*y(2)-x(3)*y(3)-x(4)*y(4); ...
   x(1)*y(2)+x(2)*y(1)+x(3)*y(4)-x(4)*y(3); ...
   x(1)*y(3)-x(2)*y(4)+x(3)*y(1)+x(4)*y(2); ...
   x(1)*y(4)+x(2)*y(3)-x(3)*y(2)+x(4)*y(1)];
end
