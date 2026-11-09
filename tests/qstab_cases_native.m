function tests = qstab_cases_native()
%QSTAB_CASES_NATIVE Integration tests that require MathWorks quaternion.
    tests={...
        'Native quaternion shape and coefficient components',@shape;...
        'Native scalar algebraic roots and norm polynomial',@polynomial;...
        'Native right eigenvalues of a noncommuting triangular matrix',@matrix;...
        'Native componentwise interval robust stability',@interval;...
        'Native left/right critical-parameter frequency witnesses',@sides;...
        'Native delay coefficients and declared grade',@delay;...
        'Native quaternion polytope uses the four-generator dispatcher',@fourGenerator;...
        'Native delay four-generator witness preserves input coefficients',@fourGeneratorDelay;...
        'Native varying-leading polytope has a full-hull certificate',@varyingLeading;...
        'Native discrete coefficient intervals retain shared coefficient parameters',@discreteInterval;...
        'Native Kharitonov cube and paired faces use ordinary public family data',@cubeReduction;...
        'Native quadratic box detects its unstable two-dimensional interior',@quadraticInterior};
end
function shape()
    a=reshape([2 3 1],1,1,3); b=zeros(size(a)); c=b; d=b; c(1)=1;
    P=quaternion(a,b,c,d); C=qstab.polycomponents(P);
    assert(isequal(size(P),[1 1 3]));
    assert(isequal(C,[2 3 1;0 0 0;1 0 0;0 0 0]));
    [aa,bb,cc,dd]=qstab.components(P);
    assert(isequal(aa,a)&&isequal(bb,b)&&isequal(cc,c)&&isequal(dd,d));
end
function polynomial()
    P=quaternion(reshape([0 1],1,1,2),zeros(1,1,2),reshape([-1 0],1,1,2),zeros(1,1,2));
    [z,I]=qstab.roots(P); assert(I.isValid&&abs(z-1i)<1e-12);
    assert(isequal(qstab.normpoly(P),[1 0 1]));
    j=quaternion(0,0,1,0); i=quaternion(0,1,0,0);
    [a,b,c,d]=parts(i-j); assert(norm([a b c d])>1);
end
function matrix()
    A=quaternion([-1 0;0 -2],[1 1;0 0],[0 2;0 1],[0 3;0 0]);
    [z,I]=qstab.eig(A); assert(I.isValid);
    assert(max(abs(sort(real(z))-[-2;-1]))<1e-12&&max(abs(imag(z)-1))<1e-12);
    [~,R] = qstab.isstable(A,'Type','matrix','Domain','continuous');
    assert(strcmp(R.status,'stable')&&~R.spectrumComputed);
end
function interval()
    L=quaternion([1 .999],[-.02 0],[-.02 0],[-.02 0]);
    U=quaternion([1.1 1.001],[.02 0],[.02 0],[.02 0]);
    [~,R] = qstab.robuststable({L,U},'Uncertainty','interval','Variable','s');
    assert(strcmp(R.status,'stable')&&R.certificateComplete);
end
function sides()
    p0=quaternion([1 2 1],[0 0 0],[1 0 0],[0 0 0]);
    p1=quaternion([0 0],[0 1],[0 0],[1 0]);
    [rL,IL]=qstab.freqmargin(p0,p1,0,'Side','left');
    [rR,IR]=qstab.freqmargin(p0,p1,sqrt(3)-1,'Side','right');
    assert(IL.valid&&IR.valid&&abs(rL-sqrt(2))<1e-12);
    assert(abs(rR-sqrt(2)*(2-sqrt(3)))<1e-12);
end
function delay()
    P=quaternion([1 -.5 0],[0 0 0],[0 0 0],[0 0 0]);
    [~,R] = qstab.isstable(P,'Variable','zinv'); assert(strcmp(R.status,'stable'));
    [z,I]=qstab.roots(P,'Variable','zinv'); assert(I.isValid);
    assert(I.infiniteCount==1&&abs(z(isfinite(z))-2)<1e-12);
end

function fourGenerator()
    V=cell(1,5);
    for j=1:5
        V{j}=quaternion([1+.1*j,1],[.02*j,0],[.01*(-1)^j,0],[.02,0]);
    end
    [~,R] = qstab.robuststable(V,'Method','zeroexclusion', ...
        'Reduction','four-generator','MaxTime',Inf);
    assert(strcmp(R.status,'stable') && R.reduction.stableSubfamilies==5);
    assert(R.certificateComplete && ~R.validated);
end
function fourGeneratorDelay()
    V=repmat({quaternion([1,2,0],[0,0,0],[0,0,0],[0,0,0])},1,5);
    [~,R] = qstab.zeroexclude(V,'Variable','zinv','Reduction','four-generator','MaxTime',Inf);
    assert(strcmp(R.status,'notStable') && numel(R.witness.originalWeights)==5);
    assert(isequal(R.witness.originalComponents(1,:),[1,2,0]));
    assert(isequal(R.witness.components(1,:),[0,2,1]));
end

function varyingLeading()
    V={quaternion([1,1],[1,1],[0,0],[0,0]),quaternion([2,2],[-1,-1],[0,0],[0,0])};
    [~,R] = qstab.zeroexclude(V,'MaxTime',Inf);
    assert(strcmp(R.status,'stable')&&R.leading.isSeparated);
end
function discreteInterval()
    L=quaternion([-.3,1],[-.05,0],[-.1,0],[-.05,0]);
    U=quaternion([-.1,1],[.05,0],[.1,0],[.05,0]);
    [~,R] = qstab.zeroexclude({L,U},'Uncertainty','interval','Variable','z','MaxTime',Inf);
    assert(strcmp(R.status,'stable')&&R.certificateComplete);
end
function cubeReduction()
    L=quaternion([1,1],[0,0],[0,0],[0,0]);
    U=quaternion([2,2],[0,0],[0,0],[0,0]);
    [K,I]=qstab.kharitonovcube({L,U});assert(isstruct(K)&&I.activeParameterCount==2);
    [~,R] = qstab.zeroexclude({L,U},'Uncertainty','interval','Reduction','kharitonov-faces','MaxTime',Inf);
    assert(strcmp(R.status,'stable')&&R.reduction.totalSubfamilies==4);
end
function quadraticInterior()
    L=quaternion([0,1,1],[2,1,0],[0,-2,0],[0,-2,0]);
    U=quaternion([0,1,1],[2,1,0],[0,2,0],[0,2,0]);
    [~,R] = qstab.zeroexclude({L,U},'Uncertainty','interval','Reduction','quadratic-faces','MaxTime',Inf);
    assert(strcmp(R.status,'notStable')&&~R.isStable);
end
