function tests=qstab_cases_oneparameter_numeric()
%QSTAB_CASES_ONEPARAMETER_NUMERIC Ordered witnesses, geometry and time units.
tests={ ...
 'E1 ordered Hamilton and SVD margins retain both evaluation witnesses',@ordered; ...
 'E2 rank-one perturbation needs no inverse and has a unique gain',@rankOne; ...
 'E3 nearest-gain sphere and repeated-value singleton are distinguished',@geometry; ...
 'Rotated coefficient orientations preserve the complete gain sphere',@rotatedGeometry; ...
 'E4 near multiplicity is a poorly conditioned point, not a sphere',@nearRepeated; ...
 'E5 tangent candidates remain unresolved in the numerical ray method',@tangencies; ...
 'E6 flat frequency minimum retains a valid witness',@flat; ...
 'E7 frequency scaling preserves original frequencies, gains and ray parameters',@units; ...
 'E7 global numerical brackets and frequency cells retain original units',@globalUnits; ...
 'A nonmonic quaternion factor preserves the ordered parameter transformation',@nonmonic; ...
 'Family entry points forward scaling and preserve the caller inputs',@familyEntry; ...
 'Scaling options reject incompatible domains and unrepresentable coefficients',@invalid ...
};
end

function ordered()
[p0,p1]=example();w=[0;sqrt(3)-1;1];
for side={'left','right'}
    [r,I]=qstab.freqmargin(p0,p1,w,'Side',side{1});
    assert(isequal(size(r),size(w))&&all(I.valid));
    if strcmp(side{1},'left'),target=sqrt(2+w.^4./(1+w).^2);
    else,target=sqrt((w.^4+2*w.^2-4*w+2)./(1+w).^2);end
    assert(max(abs(r-target))<3e-12&&max(I.svdRelativeDifference)<1e-12);
    assert(strcmp(I.geometry{1}.kind,'point')&&I.frequencyScale==1);
    assert(~I.scaling.applied);
    for k=1:numel(w)
        Z=I.zeroWitnesses{k};assert(strcmp(Z.status,'computed'));
        assert(max([Z.rightResidual Z.leftResidual Z.boundaryClassResidual])<1e-11);
        assert(abs(norm(Z.rightEvaluationComponents)-w(k))<1e-11);
        assert(abs(norm(Z.leftEvaluationComponents)-w(k))<1e-11);
    end
end
end

function rankOne()
[r,I]=qstab.freqmargin([1.5 1.5 1],[1i 1],1);
assert(I.valid&&abs(r-sqrt(5/8))<1e-12);
assert(norm(I.parameterComponents-[-.75;.25;0;0])<1e-12);
assert(strcmp(I.geometry{1}.kind,'point')&&I.singularValueGapRatio>.99);
assert(norm(I.zeroWitnesses{1}.rightEvaluationComponents-[0;1;0;0])<1e-12);
end

function geometry()
[r,I]=qstab.freqmargin([1 1 1],1,1/sqrt(2));G=I.geometry{1};
assert(abs(r-sqrt(3)/2)<1e-12&&strcmp(G.kind,'sphere'));
assert(norm(G.centerComponents-[-.5;0;0;0])<1e-12);
assert(abs(G.sphereRadius-1/sqrt(2))<1e-12&&G.dimension==2);
assert(norm(G.basisComponents'*G.basisComponents-eye(3),'fro')<1e-12);
assert(G.boundaryResidual<1e-12&&G.selectedMembershipResidual<1e-12);
[r,J]=qstab.freqmargin([1 2 1],[0 1],1);G=J.geometry{1};
assert(abs(r-2)<1e-12&&strcmp(G.kind,'point'));
assert(norm(G.centerComponents-[-2;0;0;0])<1e-12&&G.sphereRadius==0);
[~,J]=qstab.freqmargin([1 1 1],1,0);
assert(strcmp(J.geometry{1}.kind,'point'));
[r,J]=qstab.freqmargin([1 3 3 1],[1 0 1],1);
assert(isinf(r)&&J.valid&&strcmp(J.geometry{1}.kind,'empty'));
end

function rotatedGeometry()
q0=[1;2;-1;3];q1=[2;-1;1;1];
C0=q0*[1 1 1];C1=q1;
for side={'left','right'}
    [r,I]=qstab.freqmargin(parts(C0),parts(C1),1/sqrt(2),'Side',side{1});G=I.geometry{1};
    assert(I.valid&&strcmp(G.kind,'sphere'));
    assert(abs(r-sqrt(3)/2*norm(q0)/norm(q1))<3e-12);
    assert(abs(G.sphereRadius-norm(q0)/norm(q1)/sqrt(2))<3e-12);
    assert(G.boundaryResidual<1e-12&&G.selectedMembershipResidual<1e-12);
end
end

function nearRepeated()
e=1e-10;w=1/sqrt(2)+e/2;
C=[1 1 1;0 0 0;e 0 0;0 0 0];
[r,I]=qstab.freqmargin(parts(C),1,w);G=I.geometry{1};
assert(I.valid&&strcmp(G.kind,'point')&&G.inverseGapEstimate>1e8);
assert(norm(I.parameterComponents-[w*w-1;0;w-e;0])<2e-12);
assert(abs(r-sqrt((1-w*w)^2+(w-e)^2))<1e-12);
% A nonzero perturbation far below the numerical gap resolution is NOT
% reclassified as an exact repeated-value sphere.
C(3,1)=1e-17;[~,J]=qstab.freqmargin(parts(C),1,1/sqrt(2));
assert(strcmp(J.geometry{1}.kind,'undetermined'));
end

function tangencies()
[p0,~]=example();
for C1={[-1 -1;0 0;-1 -1;0 0],[-1 -1;0 1;-1 -1;0 0]}
    [ab,I]=qstab.raymargin(p0,parts(C1{1}),1);
    assert(any(strcmp(I.status,{'undetermined','partial'}))&&isnan(ab(2)));
    assert(abs(I.candidateInterval(2)-1)<1e-5);
end
end

function flat()
e=1e-8;w=sqrt(e);[r,I]=qstab.freqmargin([2+e 2 1],1,w);
assert(I.valid&&abs(r-sqrt(4+4*e))<2e-12&&I.svdRelativeDifference<1e-12);
end

function units()
[p0,p1]=example();[K0,K1]=scaledExample();w=sqrt(3)-1;u=parts([0;1;0;1]);
for side={'left','right'}
    [r,I]=qstab.freqmargin(p0,p1,w,'Side',side{1});
    [rk,K]=qstab.freqmargin(K0,K1,1e6*w,'Side',side{1});
    assert(K.valid&&K.scaling.applied&&K.frequencyScale>1e5);
    assert(abs(r-rk)<3e-12&&norm(I.parameterComponents-K.parameterComponents)<3e-12);
    assert(norm(K.boundaryComponents/1e6-I.boundaryComponents)<3e-12);
    Z=K.zeroWitnesses{1};assert(strcmp(Z.status,'computed'));
    assert(abs(norm(Z.rightEvaluationComponents)/1e6-w)<3e-12);
    [rk,J]=qstab.freqmargin(K0,K1,1e6*w,'Side',side{1},'FrequencyScale',1e6);
    assert(J.valid&&abs(rk-r)<3e-12&&abs(J.workingFrequency-w)<eps);
    [ab,A]=qstab.raymargin(p0,p1,u,'Side',side{1});
    [ak,B]=qstab.raymargin(K0,K1,u,'Side',side{1});
    assert(strcmp(B.status,'computed')&&abs(ak(2)-ab(2))<2e-11);
    assert(isequal(A.directionComponents,B.directionComponents)&&B.scaling.parameterUnchanged);
end
end

function globalUnits()
[p0,p1]=example();[K0,K1]=scaledExample();
opts={'Side','right','MaxCells',1,'MaxTime',Inf,'AbsTol',1e-8,'RelTol',1e-8};
[ab,A]=qstab.stabilityradius(p0,p1,opts{:},'FrequencyScale',1);
[ak,K]=qstab.stabilityradius(K0,K1,opts{:},'FrequencyScale',1e6);
assert(max(abs(ab-ak))<2e-11&&~K.validated);
assert(abs(K.frequency/1e6-A.frequency)<2e-11);
assert(norm(K.boundaryComponents/1e6-A.boundaryComponents)<2e-11);
assert(max(max(abs(K.diagnostics.cells(:,1:2)/1e6-A.diagnostics.cells(:,1:2))))<2e-11);
assert(ak(1)<=sqrt(2)*(2-sqrt(3))&&ak(2)>=sqrt(2)*(2-sqrt(3)));
end

function nonmonic()
[p0,p1]=example();C0=qstab.polycomponents(p0);C1=qstab.polycomponents(p1);
h=[1;2;-1;.5];hi=[h(1);-h(2:4)]/(h'*h);u=[0;1;0;1];
D0=leftmultiply(h,C0);D1=leftmultiply(h,C1);w=sqrt(3)-1;
for side={'left','right'}
    [r,I]=qstab.freqmargin(p0,p1,w,'Side',side{1});
    [rr,J]=qstab.freqmargin(parts(D0),parts(D1),w,'Side',side{1});
    if strcmp(side{1},'left'),expected=mul(mul(h,I.parameterComponents),hi);d=mul(mul(h,u),hi);
    else,expected=I.parameterComponents;d=u;end
    assert(J.valid&&abs(rr-r)<3e-12&&norm(J.parameterComponents-expected)<3e-12);
    [ab,~]=qstab.raymargin(p0,p1,parts(u),'Side',side{1});
    [bb,R]=qstab.raymargin(parts(D0),parts(D1),parts(d),'Side',side{1});
    assert(strcmp(R.status,'computed')&&abs(ab(2)-bb(2))<3e-11);
end
end

function familyEntry()
[p0,p1]=scaledExample();F=qstab.family('oneparameter',p0,p1,'Side','right');before=F;
[r,I]=qstab.freqmargin(F,1e6*(sqrt(3)-1));
assert(I.valid&&abs(r-sqrt(2)*(2-sqrt(3)))<1e-12&&isequaln(F,before));
assert(strcmp(I.parameterDomain.side,'right'));
[a,A]=qstab.raymargin(F,parts([0;1;0;1]));
assert(strcmp(A.status,'computed')&&abs(a(2)-(2-sqrt(3)))<1e-11);
end

function invalid()
bad(@()qstab.freqmargin([1 1],1,0,'FrequencyScale',0),'qstab:oneparameter:FrequencyScale');
bad(@()qstab.freqmargin([1 1],1,0,'FrequencyScale','bad'),'qstab:oneparameter:FrequencyScale');
bad(@()qstab.freqmargin([-.5 1],1,0,'Variable','z','FrequencyScale',2),'qstab:freqmargin:FrequencyScale');
bad(@()qstab.raymargin(-1,1,1,'Type','matrix','FrequencyScale',2),'qstab:raymargin:FrequencyScale');
bad(@()qstab.stabilityradius([-.5 1],1,'Variable','z','FrequencyScale',2),'qstab:stabilityradius:FrequencyScale');
bad(@()qstab.freqmargin([1 2 1],1,0,'FrequencyScale',1e-310),'qstab:oneparameter:ScalingRange');
end

function bad(f,id)
try,f();catch e,assert(strcmp(e.identifier,id));return;end
error('qstab:test:ExpectedError','Expected %s.',id);
end
function [p0,p1]=example()
p0=parts([1 2 1;0 0 0;1 0 0;0 0 0]);p1=parts([0 0;0 1;0 0;1 0]);
end
function [p0,p1]=scaledExample()
p0=parts([1e80 2e74 1e68;0 0 0;1e80 0 0;0 0 0]);p1=parts([0 0;0 1e74;0 0;1e80 0]);
end
function p=parts(C)
p=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function C=leftmultiply(q,C)
for j=1:size(C,2),C(:,j)=mul(q,C(:,j));end
end
function q=mul(a,b)
q=[a(1)*b(1)-a(2:4)'*b(2:4);a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
end
