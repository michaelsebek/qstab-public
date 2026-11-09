function cases = qstab_cases_leading()
%QSTAB_CASES_LEADING Full-family degree certificates, not endpoint heuristics.
cases={ ...
 'Common nonreal leading coefficient has its quaternion norm bound',@common; ...
 'Cartesian leading boxes have exact distance and exact zero witnesses',@boxes; ...
 'A separating plane succeeds although the coordinate box contains zero',@plane; ...
 'Correlated affine leading coefficients retain their common parameter',@correlation; ...
 'Five-generator interior degree loss is caught before four-generator tests',@five; ...
 'A tiny nonzero leading distance is not misreported as a zero witness',@nearZero; ...
 'Varying-leading Hurwitz polytopes are tested without memberwise normalization',@varyingCT; ...
 'Varying-leading Schur and delay polytopes retain equivalent coefficients',@varyingDT; ...
 'Four-generator reduction accepts a certified varying-leading full hull',@four; ...
 'Leading iteration budget yields undetermined, not a false certificate',@budget; ...
 'Leading-bound options and invalid inputs are rejected explicitly',@invalid; ...
 'Leading bounds scale with a common positive coefficient factor',@scaling};
end
function common()
P=poly([1,1;2,2;3,3;4,4]);[~,R] = qstab.leadingbound({P,P});
assert(R.isSeparated&&abs(R.lowerBound-sqrt(30))<1e-13&&~R.validated);
end
function boxes()
L=poly([1,1;0,-2;0,0;0,0]);U=poly([1,2;0,-1;0,0;0,0]);
[~,R] = qstab.leadingbound({L,U},'Uncertainty','interval');
assert(R.isSeparated&&abs(R.lowerBound-sqrt(2))<1e-13);
[~,R] = qstab.leadingbound({[1,-1],[1,1]},'Uncertainty','interval');
assert(strcmp(R.status,'zeroWitness')&&all(R.witness.components(:,end)==0));
end
function plane()
V={poly([1,1;0,0;0,0;0,0]),poly([0,0;1,1;0,0;0,0])};
[~,R] = qstab.leadingbound(V);
assert(R.isSeparated&&R.lowerBound>0.7&&R.lowerBound<=1/sqrt(2));
end
function correlation()
F=struct('center',poly([1,1;1,1;0,0;0,0]), ...
 'generators',{{poly([1,1;-1,-1;0,0;0,0])}},'bounds',[-2,2]);
[~,R] = qstab.leadingbound(F);
assert(R.isSeparated&&R.lowerBound>1.4);
end
function five()
Q=[eye(4),-ones(4,1)];V=cell(1,5);
for j=1:5,V{j}=poly([Q(:,j),Q(:,j)]);end
[~,L] = qstab.leadingbound(V);assert(strcmp(L.status,'zeroWitness'));
[~,R] = qstab.zeroexclude(V,'Reduction','four-generator','MaxTime',Inf);
assert(strcmp(R.status,'assumptionFailed')&&R.reduction.processedSubfamilies==0);
assert(numel(R.witness.originalWeights)==5&&all(R.witness.originalComponents(:,end)==0));
end
function nearZero()
V={poly([1e-13,1e-13;1,1;0,0;0,0]),poly([1e-13,1e-13;-1,-1;0,0;0,0])};
[~,R] = qstab.leadingbound(V);assert(strcmp(R.status,'undetermined')&&isempty(R.witness));
end
function varyingCT()
V={poly([1,1;1,1;0,0;0,0]),poly([2,2;-1,-1;0,0;0,0])};
[~,R] = qstab.zeroexclude(V,'MaxTime',Inf);stable(R);
assert(R.leading.isSeparated&&R.leading.lowerBound>0);
end
function varyingDT()
V={poly([-.2,1;-.2,1;0,0;0,0]),poly([-.4,2;.2,-1;0,0;0,0])};
[~,R] = qstab.zeroexclude(V,'Variable','z','MaxTime',Inf);stable(R);
W=cell(size(V));for j=1:numel(V),C=qstab.polycomponents(V{j});W{j}=poly(C(:,end:-1:1));end
[~,S] = qstab.zeroexclude(W,'Variable','zinv','MaxTime',Inf);stable(S);
assert(R.processedCells==S.processedCells&&R.frequencyLimit==S.frequencyLimit);
end
function four()
V=cell(1,5);for j=1:5,V{j}=poly([1+j,1+j;.1*j,.1*j;0,0;0,0]);end
[~,R] = qstab.zeroexclude(V,'Reduction','four-generator','MaxTime',Inf);stable(R);
assert(R.reduction.stableSubfamilies==5&&R.leading.isSeparated);
end
function budget()
V={poly([1,1;0,0;0,0;0,0]),poly([-2,-2;1,1;0,0;0,0])};
[~,A] = qstab.leadingbound(V,'MaxIterations',0);assert(strcmp(A.status,'undetermined'));
[~,B] = qstab.leadingbound(V,'MaxIterations',256);assert(B.isSeparated&&B.lowerBound>0);
end
function invalid()
bad={-1,NaN,1.5,Inf};
for j=1:numel(bad)
 caught=false;try,qstab.leadingbound({[1,1]},'MaxIterations',bad{j});catch,caught=true;end
 assert(caught);
end
caught=false;try,qstab.leadingbound({[1,1]},'Variable','unknown');catch,caught=true;end
assert(caught);
end
function scaling()
[~,A] = qstab.leadingbound({[1,1],[2,2]});[~,B] = qstab.leadingbound({[1e4,1e4],[2e4,2e4]});
assert(A.isSeparated&&B.isSeparated&&abs(B.lowerBound/1e4-A.lowerBound)<1e-12);
end
function P=poly(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
 'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function stable(R)
assert(strcmp(R.status,'stable')&&R.certificateComplete&&~R.validated,R.message);
end
