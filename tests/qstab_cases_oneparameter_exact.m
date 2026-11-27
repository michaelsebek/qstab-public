function cases=qstab_cases_oneparameter_exact()
%QSTAB_CASES_ONEPARAMETER_EXACT Optional rational and algebraic proofs.
% The harness should skip this group when Symbolic Math Toolbox is absent.
% A successful run recomputes proofs; stored paper results are not evidence.
cases={ ...
 'Exact radius bounds cover all frequencies and reject false certificates',@radiusBasic; ...
 'Exact left/right radii preserve noncommutative multiplication order',@radiusSides; ...
 'Tangencies and norm-induced multiplicity retain the first exact contact',@tangencies; ...
 'An algebraic ray contact lifts both evaluation zeros exactly',@algebraicLifting; ...
 'No-positive-root and zero-perturbation proofs return infinite endpoints',@infinityCases; ...
 'A later real-pair guardian contact cannot replace the first contact',@laterContact; ...
 'Directions are not normalized and precision budgets remain explicit',@directionAndBudget; ...
 'Binary rationals, explicit decimals and trailing-zero private padding are distinct',@rationalInputs; ...
 'Exact certificates check nominal stability, degree and family-domain assumptions',@assumptions; ...
 'Exact reports remain ordinary saveable MATLAB data',@ordinaryReports; ...
 'Several simultaneous imaginary classes have an exact implicit lifting',@multipleClasses; ...
 'Decimal and scientific text denotes exact fractions without binary rounding',@exactText; ...
 'Approximate symbolic numbers and malformed rational text are rejected',@invalidExactData; ...
 'Decimal and fraction radius bounds give identical exact Sturm certificates',@decimalRadiusEquivalence; ...
 'Polynomial identities and certificates are quiet and independent of symbolic precision',@quietExactCertificates};
end

function radiusBasic()
[b,I]=qstab.radiuscertificate([1 1 1],1,{'0.8660','0.8661'});
assert(I.certified&&I.ok&&I.validated&&b(1)<sqrt(3)/2&&b(2)>sqrt(3)/2, ...
    'The valid decimal bracket [0.8660,0.8661] was not certified (status: %s).',I.status);
assert(I.lowerProof.nonnegativeRootCount==0&&I.upperProof.nonnegativeRootCount>0, ...
    'The valid bracket must have no lower-bound root and a positive upper-bound root count.');
[b,I]=qstab.radiuscertificate([1 1 1],1,{'0.9','1'});
assert(~I.certified&&all(isnan(b))&&~I.lowerProof.accepted, ...
    'The false lower bound 0.9 for sqrt(3)/2 must be rejected.');
[b,I]=qstab.radiuscertificate([1 1],1,{'0.5','0.75'});
assert(~I.certified&&all(isnan(b))&&~I.upperProof.accepted, ...
    'The false upper bound 0.75 for radius 1 must be rejected.');
[b,I]=qstab.radiuscertificate([1 1],1,{'0.5','1'},'UpperX','0');
assert(I.certified&&b(2)==1&&strcmp(I.upperProof.Qvalue,'0'), ...
    'The exact zero-frequency witness for radius 1 must be retained.');
end
function radiusSides()
[P0,P1]=separating();
[b,L]=qstab.radiuscertificate(P0,P1,{'1.414213562373','1.414213562374'},'Side','left','UpperX','0');
assert(L.certified&&b(1)<=sqrt(2)&&b(2)>=sqrt(2));
[b,R]=qstab.radiuscertificate(P0,P1,{'0.378937381963','0.378937381964'},'Side','right');
assert(R.certified&&b(1)<=sqrt(2)*(2-sqrt(3))&&b(2)>=sqrt(2)*(2-sqrt(3)));
end
function tangencies()
[P0,~]=separating();
P1=parts([-1 -1;0 0;-1 -1;0 0]);
[b,I]=qstab.raycertificate(P0,P1,1);
assert(I.certified&&isequal(b,[1 1])&&I.lifting.zeroContact&&I.lifting.exact);
P1=parts([-1 -1;0 1;-1 -1;0 0]);
[b,I]=qstab.raycertificate(P0,P1,1);
assert(I.certified&&isequal(b,[1 1])&&I.lifting.zeroContact);
[b,I]=qstab.raycertificate([3 2 2 1],[6 1 1],1);
assert(I.certified&&isequal(b,[1 1])&&I.lifting.exact&&I.lifting.sphereOfEvaluationZeros);
end
function algebraicLifting()
[P0,P1]=separating();d=parts([0;1;0;1]);
[b,I]=qstab.raycertificate(P0,P1,d,'Side','right','AbsTol',1e-8);
assert(I.certified&&b(1)<=2-sqrt(3)&&b(2)>=2-sqrt(3)&&I.widthTargetMet);
assert(I.lifting.exact&&I.lifting.rightResidualExactlyZero&&I.lifting.leftResidualExactlyZero);
lo=sym(I.exactBounds{1});hi=sym(I.exactBounds{2});
assert(isAlways(sym(b(1),'f')<=lo)&&isAlways(sym(b(2),'f')>=hi));
end
function infinityCases()
[b,I]=qstab.raycertificate([1 1],1,1);assert(I.certified&&all(isinf(b)));
[b,I]=qstab.raycertificate([1 1],0,1);assert(I.certified&&all(isinf(b)));
[b,I]=qstab.radiuscertificate([1 1],0,[0 1]);assert(I.certified&&all(isinf(b)));
[P0,P1]=separating();[b,I]=qstab.raycertificate(P0,P1,1);
assert(I.certified&&all(isinf(b))&&I.isolation.positiveRootCount==0);
end
function laterContact()
[b,I]=qstab.raycertificate([2 3 1],[-1 -1],1);
assert(I.certified&&isequal(b,[2 2])&&I.lifting.zeroContact);
assert(I.isolation.positiveRootCount==2);
end
function directionAndBudget()
[b,I]=qstab.raycertificate([1 1],-1,2);assert(I.certified&&isequal(b,[.5 .5]));
assert(strcmp(I.directionNormSquared,'4'));
[P0,P1]=separating();d=parts([0;1;0;1]);
[b,I]=qstab.raycertificate(P0,P1,d,'Side','right','MaxBisections',1,'VerifyLifting',false);
assert(~I.widthTargetMet);
if I.certified,assert(b(1)<=2-sqrt(3)&&b(2)>=2-sqrt(3));else,assert(all(isnan(b)));end
end
function rationalInputs()
C0=[.1 1;0 0;0 0;0 0];C1=[-1;0;0;0];
[~,B]=qstab.raycertificate([.1 1],-1,1,'VerifyLifting',false);
assert(strcmp(B.inputSemantics,'exact-binary-input'));
assert(~strcmp(B.inputComponents.P0{1},'1/10'));
R=struct('P0',{{'1/10','1';'0','0';'0','0';'0','0'}},'P1',C1);
[b,I]=qstab.raycertificate([.1 1],[-1 0],1,'RationalData',R,'VerifyLifting',false);
assert(I.certified&&strcmp(I.inputComponents.P0{1},'1/10')&&b(1)<=.1&&b(2)>=.1);
R.P0=C0;R.P0(1,1)=.2;mustError(@()qstab.raycertificate([.1 1],-1,1,'RationalData',R),'qstab:exact:RationalMismatch');
R.P0={'1/10','1';'1e-400','0';'0','0';'0','0'};
mustError(@()qstab.raycertificate([.1 1],-1,1,'RationalData',R),'qstab:exact:RationalMismatch');
end
function assumptions()
[b,I]=qstab.raycertificate([-1 1],1,1);assert(~I.certified&&strcmp(I.status,'assumptionFailed')&&all(isnan(b)));
[b,I]=qstab.radiuscertificate([1 1],[1 1],[0 1]);assert(~I.certified&&strcmp(I.status,'assumptionFailed')&&all(isnan(b)));
mustError(@()qstab.raycertificate([1 1],1,0),'qstab:exact:Direction');
mustError(@()qstab.raycertificate([1 1],1,1,'Variable','z'),'qstab:exact:Variable');
F=qstab.family('oneparameter',[1 1],-1,'Bounds',[-1 1;0 0;0 0;0 0]);
mustError(@()qstab.raycertificate(F,1),'qstab:family:BoundedMargin');
[b,I]=qstab.raycertificate(F,1,'ExtendParameter',true);
assert(I.certified&&isequal(b,[1 1])&&I.parameterDomain.explicitExtension);
end
function ordinaryReports()
[~,I]=qstab.radiuscertificate([1 1],1,{'0.9','1.1'});
assert(~hasSymbolic(I));
[~,J]=qstab.raycertificate([1 1],-1,1);assert(~hasSymbolic(J));
end
function multipleClasses()
% At gamma=1: (s^2+1)*(s^2+4), with two distinct imaginary classes.
[b,I]=qstab.raycertificate([4 3 5 1 1],[0 -3 0 -1],1);
assert(I.certified&&isequal(b,[1 1])&&I.lifting.exact);
assert(strcmp(I.lifting.status,'verifiedImplicit')&&I.lifting.positiveFrequencyCount==2);
assert(I.lifting.normIdentitiesVerified&&~I.lifting.numericalWitness);
end
function exactText()
% The expected values are built from integers, independently of the text parser.
input={'0.8660','.125e1','-1.25e+2 / +2.5e-1'; ...
       '+00012.3400e-2','9007199254740993.125','-0.000e+400'};
expected=[sym(433)/500,sym(5)/4,sym(-500); ...
          sym(617)/5000,sym('72057594037927945')/8,sym(0)];
actual=qstab.internal.exactOneParameter('rational',input);
assert(isequal(actual,expected), ...
    'Decimal/scientific/fraction text must preserve its exact rational value and array shape.');
tiny=qstab.internal.exactOneParameter('rational','1e-400');
assert(isequal(tiny*sym(10)^400,sym(1)), ...
    'A nonzero exact decimal below double range must not underflow to zero.');
large=qstab.internal.exactOneParameter('rational','1e400');
assert(isequal(large,sym(10)^400), ...
    'An exact decimal above double range must remain a finite symbolic rational.');
decimal=qstab.internal.exactOneParameter('rational','0.1');
binary=qstab.internal.exactOneParameter('rational',0.1);
assert(isequal(decimal,sym(1)/10)&&~isequal(decimal,binary), ...
    'Decimal text 0.1 and the stored binary double 0.1 have different exact semantics.');
end
function invalidExactData()
for value={'1/0','1/-0.00e-30','NaN','Inf','sqrt(2)','1+2','3/4/5','1e'}
    mustError(@()qstab.internal.exactOneParameter('rational',value{1}),'qstab:exact:Rational');
end
mustError(@()qstab.internal.exactOneParameter('rational',vpa(sym(1)/3,20)), ...
    'qstab:exact:Rational');
mustError(@()qstab.internal.exactOneParameter('rational',vpa(sym(3),20)), ...
    'qstab:exact:Rational');
mustError(@()qstab.internal.exactOneParameter('rational',sqrt(sym(2))), ...
    'qstab:exact:Rational');
mustError(@()qstab.internal.exactOneParameter('rational',sym('qstabExactTestFreeVariable')), ...
    'qstab:exact:Rational');
% The public route must enforce the same contract as the internal parser.
mustError(@()qstab.radiuscertificate([1 1 1],1,vpa([sym(433)/500,sym(8661)/10000],20)), ...
    'qstab:exact:Rational');
x=sym('qstabExactTestFloatPolynomial');
mustError(@()qstab.internal.exactOneParameter('sturm',vpa(sym(1)/3,20)*x+1,x), ...
    'qstab:exact:Rational');
end
function decimalRadiusEquivalence()
[bd,D]=qstab.radiuscertificate([1 1 1],1,{'0.8660','8.661e-1'});
[bf,F]=qstab.radiuscertificate([1 1 1],1,{'433/500','8661/10000'});
assert(D.certified&&F.certified&&isequal(bd,bf), ...
    'Equivalent decimal and fraction brackets must both certify sqrt(3)/2.');
assert(isequal(D.exactBounds,{'433/500','8661/10000'})&&isequal(D.exactBounds,F.exactBounds), ...
    'The certificate must record reduced integer fractions, not symbolic floating-point endpoints.');
assert(isequal(D.polynomials,F.polynomials)&& ...
       isequal(D.lowerProof,F.lowerProof)&&isequal(D.upperProof,F.upperProof), ...
    'Equivalent exact inputs must produce the same polynomial and Sturm evidence.');
assert(D.lowerProof.nonnegativeRootCount==0&&D.upperProof.nonnegativeRootCount==2, ...
    'Square-free Sturm counts must retain both positive contacts of the upper radius bound.');
end
function quietExactCertificates()
% A warning promoted to an error catches the native MATLAB regression.
id='symbolic:sym:isAlways:TruthUnknown';priorWarning=warning('query',id);
priorDigits=digits;[priorMessage,priorId]=lastwarn;
cleanup=onCleanup(@()restoreSymbolicTestState(priorWarning,priorDigits,priorMessage,priorId)); %#ok<NASGU>
warning('error',id);digits(16);
x=sym('qstabExactTestIdentity');
assert(qstab.internal.exactOneParameter('iszero',(x+1)^2-x^2-2*x-1), ...
    'An identically zero polynomial must be recognized algebraically.');
assert(~qstab.internal.exactOneParameter('iszero',x*(x-1)), ...
    'A nonzero polynomial is not an identity merely because it has individual zeros.');
[~,A]=qstab.radiuscertificate([1 1 1],1,{'0.8660','0.8661'});
digits(50);
[~,B]=qstab.radiuscertificate([1 1 1],1,{'0.8660','0.8661'});
[b,R]=qstab.raycertificate([3 2 2 1],[6 1 1],1);
assert(A.certified&&B.certified&&R.certified&&isequal(b,[1 1])&&R.lifting.exact, ...
    'The exact radius and tangency proofs must not emit TruthUnknown warnings.');
assert(isequal(A.polynomials,B.polynomials)&&isequal(A.lowerProof,B.lowerProof)&& ...
       isequal(A.upperProof,B.upperProof)&&isequal(A.exactBounds,B.exactBounds), ...
    'Exact rational certificates must be independent of the caller''s digits setting.');
currentWarning=warning('query',id);
assert(strcmp(currentWarning.state,'error')&&digits==50, ...
    'Certificate calls must preserve warning and symbolic precision settings.');
end
function restoreSymbolicTestState(w,d,m,id)
warning(w.state,w.identifier);digits(d);lastwarn(m,id);
end
function tf=hasSymbolic(x)
tf=isa(x,'sym');if tf,return;end
if isstruct(x),names=fieldnames(x);for n=1:numel(x),for k=1:numel(names),tf=tf||hasSymbolic(x(n).(names{k}));end,end
elseif iscell(x),for k=1:numel(x),tf=tf||hasSymbolic(x{k});end,end
end
function [P0,P1]=separating()
P0=parts([1 2 1;0 0 0;1 0 0;0 0 0]);P1=parts([0 0;0 1;0 0;1 0]);
end
function p=parts(C),p=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));end
function mustError(f,id)
try,f();catch e,assert(strcmp(e.identifier,id),'Unexpected error identifier: %s',e.identifier);return;end
error('qstab:test:MissingError','Expected %s.',id);
end
