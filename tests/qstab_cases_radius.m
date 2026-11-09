function tests = qstab_cases_radius()
%QSTAB_CASES_RADIUS Analytic qB01 frequency/radius regressions.
    tests={...
        'Frequency radius of s+2 with a constant quaternion perturbation',@firstOrder;...
        'Frequency sphere with identically vanishing perturbation',@emptySphere;...
        'Plane exceptional frequency of the qB01 sweep',@planeCase;...
        'Exact frequency objectives separate multiplication sides',@sideObjectives;...
        'Global first-order radius bracket and convergence',@firstOrderGlobal;...
        'Global right-side radius bracket contains the exact value',@sideGlobal;...
        'Zero perturbation and unstable-nominal assumptions',@zeroAndUnstable;...
        'A finite radius budget preserves bounds without false convergence',@budget;...
        'Common coefficient scaling preserves frequency results',@scaling;...
        'Unrepresentable radius never creates a false zero witness',@underflow;...
        'Reconstructed noncommuting witness has an algebraic boundary zero',@boundaryWitness};
end
function firstOrder()
    w=[0 .25 1 3];
    [r,I]=qstab.freqmargin([2 1],1,w);
    assert(all(I.valid)); assert(max(abs(r-sqrt(4+w.^2)))<1e-12);
    assert(max(I.evaluationResidual)<1e-12);
end
function emptySphere()
    [r,I]=qstab.freqmargin([1 3 3 1],[1 0 1],1);
    assert(isinf(r)&&I.valid); assert(I.maxPerturbationModulus==0);
end
function planeCase()
    [r,I]=qstab.freqmargin([1 2 1],[-1i 1],1);
    assert(I.valid&&abs(r-1)<1e-12);
    assert(abs(I.parameterComponents(1)+1)<1e-12);
end
function sideObjectives()
    [p0,p1]=exampleParts(); w=[0 .2 sqrt(3)-1 1 2];
    rL=qstab.freqmargin(p0,p1,w,'Side','left');
    rR=qstab.freqmargin(p0,p1,w,'Side','right');
    exactL=sqrt(2+w.^4./(1+w).^2);
    exactR=sqrt((w.^4+2*w.^2-4*w+2)./(1+w).^2);
    assert(max(abs(rL-exactL))<2e-12&&max(abs(rR-exactR))<2e-12);
end
function firstOrderGlobal()
    [~,R] = qstab.stabilityradius([1 1],1,'AbsTol',1e-3,'RelTol',1e-3,'MaxTime',Inf,'MaxCells',1024);
    assert(strcmp(R.status,'converged'));
    assert(R.lower<=1+1e-12&&R.upper>=1-1e-12&&R.gap<=.0021);
    assert(~R.validated);
    assert(R.diagnostics.bernsteinShortcut.accepted && ...
        strcmp(R.diagnostics.lowerBoundMethod,'bernstein-norm-inequality'));
    [p0,p1]=exampleParts();
    [~,L] = qstab.stabilityradius(p0,p1,'Side','left','AbsTol',5e-4,'RelTol',5e-4,...
        'MaxTime',Inf,'MaxCells',128);
    assert(strcmp(L.status,'converged') && L.diagnostics.bernsteinShortcut.accepted);
    assert(L.lower<=sqrt(2)+1e-10 && L.upper>=sqrt(2)-1e-10 && L.gap<.0013);
    assert(L.diagnostics.numberOfCells>=1 && L.diagnostics.bernsteinCells<=128);
end
function sideGlobal()
    [p0,p1]=exampleParts(); exact=sqrt(2)*(2-sqrt(3));
    [~,R] = qstab.stabilityradius(p0,p1,'Side','right','AbsTol',5e-3,'RelTol',5e-3,...
        'MaxTime',Inf,'MaxCells',512);
    assert(any(strcmp(R.status,{'converged','bounded'})));
    assert(R.lower<=exact+1e-10&&R.upper>=exact-1e-10);
    assert(R.upper-exact<2e-3);
    assert(R.diagnostics.frequencyCutoff>R.frequency);
end
function zeroAndUnstable()
    [~,R] = qstab.stabilityradius([1 1],0); assert(strcmp(R.status,'unbounded'));
    assert(isinf(R.lower)&&isinf(R.upper));
    [~,R] = qstab.stabilityradius([-1 1],1); assert(strcmp(R.status,'assumptionFailed'));
end
function budget()
    [p0,p1]=exampleParts(); exact=sqrt(2)*(2-sqrt(3));
    [~,R] = qstab.stabilityradius(p0,p1,'Side','right','AbsTol',1e-12,'RelTol',1e-12,...
        'MaxCells',1,'MaxTime',Inf);
    assert(strcmp(R.status,'bounded')); assert(R.lower<exact&&R.upper>=exact-1e-12);
    assert(R.diagnostics.numberOfCells==1);
    assert(R.diagnostics.bernsteinCells<=1 && ~R.diagnostics.bernsteinShortcut.accepted);
end
function scaling()
    for scale=[1e-200 1 1e200]
        [r,I]=qstab.freqmargin(scale*[1 1],scale,[0 1]);
        assert(all(I.valid)&&max(abs(r-[1 sqrt(2)]))<1e-12);
    end
end
function underflow()
    [r,I]=qstab.freqmargin(1e-200*[1 1],1e200,0);
    assert(~I.valid&&isnan(r));
    [~,R] = qstab.stabilityradius(1e-200*[1 1],1e200,'MaxCells',1);
    assert(strcmp(R.status,'undetermined')&&~(R.upper==0));
end
function boundaryWitness()
    [p0,p1]=exampleParts(); C0=qstab.polycomponents(p0); C1=qstab.polycomponents(p1);
    for side={'left','right'}
        w=sqrt(3)-1; [~,I]=qstab.freqmargin(p0,p1,w,'Side',side{1});
        t=I.parameterComponents; T=qstab.chi(parts(t)); C=C0;
        for k=1:size(C1,2)
            D=qstab.chi(parts(C1(:,k)));
            if strcmp(side{1},'left'), D=T*D; else, D=D*T; end
            C(:,k)=C(:,k)+[real(D(1,1));imag(D(1,1));real(D(1,2));imag(D(1,2))];
        end
        p=parts(C); [~,S] = qstab.valueset(p,1i*w);
        assert(abs(S.points)<1e-11);
    end
end
function [p0,p1]=exampleParts()
    p0=parts([1 2 1;0 0 0;1 0 0;0 0 0]);
    p1=parts([0 0;0 1;0 0;1 0]);
end
function P=parts(C)
    P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
