function tests = qstab_cases_ballfix_native()
%QSTAB_CASES_BALLFIX_NATIVE Native MathWorks quaternion regression inputs.
tests={ ...
 'Native nonreal zero-class witnesses cancel the constant exactly',@nativeZero; ...
 'Native and component CT cells retain the same exact partition',@nativePartition};
end
function nativeZero()
P=quaternion([2 .7 1],[-1 .2 0],[3 -.3 0],[.7 .1 0]);C=qstab.polycomponents(P);
for eta=[1 2 3 Inf]
    [r,I]=qstab.ballvalue(P,[.4 .2 0],0,'Exponent',eta,'Witness',true);
    W=I.witnesses{1};assert(I.valid&&W.valid&&W.evaluationResidual==0);
    assert(all(W.coefficientComponents(:,1)==0));
    assert(isequal(W.coefficientComponents(:,2:end),C(:,2:end)));
    assert(isequal(W.perturbationComponents(:,1),-C(:,1)));
    assert(abs(r-norm(C(:,1))/.4)<1e-10);
end
end
function nativePartition()
P=quaternion([2 1],[0 0],[3 0],[0 0]);Q=struct('a',[2 1],'b',[0 0],'c',[3 0],'d',[0 0]);
opts={'Exponent',3,'MaxCells',31,'MaxTime',Inf,'AbsTol',1e-12,'RelTol',1e-12};
[~,A] = qstab.ballradius(P,[1 0],opts{:});[~,B] = qstab.ballradius(Q,[1 0],opts{:});
assert(A.lower==B.lower&&A.upper==B.upper&&strcmp(A.status,B.status));
assert(isequal(A.diagnostics.cells,B.diagnostics.cells)&&A.diagnostics.partitionVerified);
T=sortrows(A.diagnostics.cells,1);
assert(T(1,1)==0&&T(end,2)==1&&all(T(1:end-1,2)==T(2:end,1)));
assert(all(T(:,2)>T(:,1))&&A.diagnostics.processedCells<=31);
end
