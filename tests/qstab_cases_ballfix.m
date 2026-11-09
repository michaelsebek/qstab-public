function tests = qstab_cases_ballfix()
%QSTAB_CASES_BALLFIX Regressions for the two 0.7.0 coefficient-ball defects.
% Existing 313 tests are retained unchanged. These checks independently test
% exact zero-frequency cancellation and a disjoint parameter partition.
tests={ ...
 'A quaternion constant cancels exactly at zero for every Holder exponent',@exactZero; ...
 'Zero-class witnesses preserve all higher coefficients under extreme scaling',@scaledZero; ...
 'A fixed nonzero constant still has infinite distance at the zero class',@fixedZero; ...
 'Requesting witnesses does not invalidate analytic nonreal class distances',@pointWithWithout; ...
 'One CT bisection produces two disjoint halves with the same split point',@ctOneSplit; ...
 'One Schur bisection produces two disjoint halves of the semicircle',@dtOneSplit; ...
 'Every bounded CT run retains an exact partition at successive budgets',@ctBudgets; ...
 'Each nonreal CT cell lower bound lies below its analytic interval minimum',@ctCellBounds; ...
 'Each nonreal Schur cell lower bound lies below its analytic interval minimum',@dtCellBounds; ...
 'A shifted CT class minimum and its witness keep the same frequency',@shiftedMinimum; ...
 'Disabling stored leaves leaves the same bounds and verified partition',@noStorage; ...
 'The preexisting Bernstein positivity shortcut retains its separate coverage',@positivityUnchanged};
end
function exactZero()
C=[2 1;0 .2;3 -.4;0 .1];P=parts(C);g=[1 .5];
for eta=[1 1.5 2 3 Inf]
    [r,I]=qstab.ballvalue(P,g,0,'Exponent',eta,'Witness',true);
    W=I.witnesses{1};assert(I.valid&&W.valid&&all(W.coefficientComponents(:,1)==0));
    assert(isequal(W.perturbationComponents(:,1),-C(:,1)));
    assert(isequal(W.coefficientComponents(:,2:end),C(:,2:end)));
    assert(W.evaluationResidual==0&&all(W.boundaryComponents==0));
    assert(abs(r-norm(C(:,1))/g(1))<1e-12&&W.fixedCoefficientsPreserved);
end
end
function scaledZero()
C=[2 .7 1;-1 .2 0;3 -.3 0;.7 .1 0];g=[.4 .2 0];
for scale=[1e-180 1 1e180]
    P=parts(scale*C);
    for eta=[1 2 Inf]
        [r,I]=qstab.ballvalue(P,scale*g,0,'Exponent',eta,'Witness',true);
        W=I.witnesses{1};assert(I.valid&&W.valid&&W.evaluationResidual==0);
        assert(all(W.coefficientComponents(:,1)==0));
        assert(isequal(W.coefficientComponents(:,2:end),scale*C(:,2:end)));
        assert(isequal(W.perturbationComponents(:,1),-scale*C(:,1)));
        assert(abs(r-norm(C(:,1))/g(1))<1e-10);
    end
end
end
function fixedZero()
P=parts([2 1;0 0;3 0;0 0]);
[r,I]=qstab.ballvalue(P,[0 1],0,'Witness',true);
assert(I.valid&&isinf(r)&&isinf(I.lower)&&isinf(I.upper)&&isempty(I.witnesses{1}));
end
function pointWithWithout()
P=parts([2 1;0 0;3 0;0 0]);z=reshape([0 1 1i 3i -1+2i 2-3i],2,3);
for eta=[1 1.5 2 3 Inf]
    [a,A]=qstab.ballvalue(P,[1 .5],z,'Exponent',eta);
    [b,B]=qstab.ballvalue(P,[1 .5],z,'Exponent',eta,'Witness',true);
    assert(all(A.valid(:))&&all(B.valid(:))&&isequal(a,b));
    if eta==1,dual=Inf;elseif isinf(eta),dual=1;else,dual=eta/(eta-1);end
    for j=1:numel(z)
        ref=hypot(real(z(j))+2,abs(imag(z(j)))-3)/norm([1 .5*abs(z(j))],dual);
        assert(abs(b(j)-ref)<1e-12&&B.witnesses{j}.valid);
    end
end
end
function ctOneSplit()
R=limited([1 1 1],[1 .2 0],3,'s');T=partition(R,1);
assert(isequal(T(:,1:2),[0 .5;.5 1]));
assert(all(T(:,3)==1)&&R.diagnostics.processedCells==3&&R.diagnostics.pointEvaluations==5);
end
function dtOneSplit()
P=parts([-.3 1;0 0;-.4 0;0 0]);R=limited(P,[1 0],3,'z');T=partition(R,pi);
assert(isequal(T(:,1:2),[0 pi/2;pi/2 pi]));
assert(all(T(:,3)==1)&&R.diagnostics.processedCells==3&&R.diagnostics.pointEvaluations==6);
end
function ctBudgets()
for budget=[1 3 7 15 31]
    R=limited([1 1 1],[1 .2 0],budget,'s');T=partition(R,1);
    assert(R.diagnostics.processedCells<=budget);
    assert(size(T,1)==(R.diagnostics.processedCells+1)/2);
    assert(R.lower<=min(T(:,4))+1e-13&&~R.validated);
end
end
function ctCellBounds()
P=parts([2 1;0 0;3 0;0 0]);R=limited(P,[1 0],63,'s');T=partition(R,1);
for j=1:size(T,1)
    lo=T(j,1)/(1-T(j,1));
    if T(j,2)==1,hi=Inf;else,hi=T(j,2)/(1-T(j,2));end
    w=min(hi,max(lo,3));ref=hypot(2,w-3);
    assert(T(j,4)<=ref+1e-10);
end
assert(R.lower<=2+1e-10&&R.upper>=2-1e-10);
end
function dtCellBounds()
a=.3;b=.4;P=parts([-a 1;0 0;-b 0;0 0]);
R=limited(P,[1 0],63,'z');T=partition(R,pi);theta=atan2(b,a);
for j=1:size(T,1)
    t=min(T(j,2),max(T(j,1),theta));ref=hypot(cos(t)-a,sin(t)-b);
    assert(T(j,4)<=ref+1e-10);
end
assert(R.lower<=.5+1e-10&&R.upper>=.5-1e-10);
end
function shiftedMinimum()
P=parts([2 1;0 0;3 0;4 0]);
[~,R] = qstab.ballradius(P,[1 0],'Exponent',3,'AbsTol',2e-3,'RelTol',2e-3, ...
    'MaxCells',2047,'MaxTime',Inf);
partition(R,1);
assert(strcmp(R.status,'converged')&&R.lower<=2+1e-10&&R.upper>=2-1e-10);
t=R.boundaryParameter;w=t/(1-t);W=R.witness;
assert(abs(w-5)<.1&&abs(R.boundaryArgument-1i*w)<1e-12);
assert(W.valid&&W.evaluationResidual<1e-8&&norm(W.boundaryComponents-[0;0;-3*w/5;-4*w/5])<1e-11);
end
function noStorage()
A=limited([1 1 1],[1 .2 0],31,'s');
[~,B] = qstab.ballradius([1 1 1],[1 .2 0],'Exponent',3,'Variable','s', ...
    'MaxCells',31,'MaxTime',Inf,'AbsTol',1e-12,'RelTol',1e-12,'StoreCells',false);
assert(strcmp(A.status,B.status)&&A.lower==B.lower&&A.upper==B.upper);
assert(B.diagnostics.partitionVerified&&B.diagnostics.completeDomainCovered);
assert(isempty(B.diagnostics.cells)&&~B.diagnostics.cellsStored);
assert(isequal(A.witness.coefficientComponents,B.witness.coefficientComponents));
end
function positivityUnchanged()
[~,R] = qstab.ballradius([1 1],[1 1],'AbsTol',2e-3,'RelTol',2e-3,'MaxCells',2047,'MaxTime',Inf);
S=R.diagnostics.positivity;assert(S.accepted&&S.attempted&&strcmp(R.status,'converged'));
T=sortrows(S.acceptedCells,1);
assert(T(1,1)==0&&T(end,2)==1&&all(T(:,2)>T(:,1)));
assert(all(T(1:end-1,2)==T(2:end,1))&&all(T(:,4)>0));
end
function R=limited(P,g,budget,variable)
[~,R] = qstab.ballradius(P,g,'Exponent',3,'Variable',variable,'MaxCells',budget, ...
    'MaxTime',Inf,'AbsTol',1e-12,'RelTol',1e-12);
end
function T=partition(R,last)
assert(R.diagnostics.partitionVerified&&R.diagnostics.completeDomainCovered);
T=sortrows(R.diagnostics.cells,1);
assert(~isempty(T)&&T(1,1)==0&&T(end,2)==last&&all(T(:,2)>T(:,1)));
assert(all(T(1:end-1,2)==T(2:end,1)));
end
function P=parts(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
