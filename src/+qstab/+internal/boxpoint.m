function E = boxpoint(M,lambda,side,tol,maxFacets,makeWitness)
%BOXPOINT Actual quaternion evaluation point, NOT a complete algebraic class.
q=double(lambda(:));n=M.grade;
E=struct('radius',NaN,'lower',0,'upper',Inf,'valid',false,'status','undetermined', ...
    'evaluationPointComponents',q,'evaluationSide',side,'gauge',[], ...
    'witness',[],'message','','powerScaleLog',0,'generatorCount',0,'validated',false);
a=norm(q);if ~isfinite(a),E.message='Quaternion argument norm overflows.';return;end
if a==0,mag=[1 zeros(1,n)];unit=[1;0;0;0];
elseif a<=1,mag=a.^(0:n);unit=q/a;
else,mag=(1/a).^(n:-1:0);unit=q/a;E.powerScaleLog=n*log(a);end
if a>0&&any(mag==0 & (any(M.Cscaled~=0,1)|any(M.Wscaled>0,1)))
    E.message='A nonzero monomial underflowed; no zero-sensitivity conclusion is made.';return
end
power=zeros(4,n+1);power(:,1)=[1;0;0;0];
for k=2:n+1,power(:,k)=localMul(power(:,k-1),unit);end
value=zeros(4,1);amplitude=0;
for k=1:n+1
    if strcmp(side,'right'),term=localMul(M.Cscaled(:,k),power(:,k));
    else,term=localMul(power(:,k),M.Cscaled(:,k));end
    value=value+mag(k)*term;amplitude=amplitude+mag(k)*norm(M.Cscaled(:,k));
end
% Collapse parallel imaginary-axis generators to the eight exact even/odd
% groups from qB02. At real points only four groups are needed.
active=find(M.Wscaled>0);map=zeros(size(M.Wscaled));sgn=zeros(size(M.Wscaled));
if a==0||all(q(2:4)==0)
    G=zeros(4,4);
    for h=1:numel(active)
        [nu,k]=ind2sub(size(M.Wscaled),active(h));
        if mag(k)==0,continue;end
        map(nu,k)=nu;sgn(nu,k)=power(1,k);
        G(nu,nu)=G(nu,nu)+M.Wscaled(nu,k)*mag(k);
    end
elseif q(1)==0
    G=zeros(4,8);
    for h=1:numel(active)
        [nu,k]=ind2sub(size(M.Wscaled),active(h));parity=mod(k-1,2);j=nu+4*parity;
        map(nu,k)=j;sgn(nu,k)=(-1)^floor((k-1)/2);
        base=zeros(4,1);base(nu)=1;
        if parity==1
            if strcmp(side,'right'),base=localMul(base,unit);else,base=localMul(unit,base);end
        end
        G(:,j)=G(:,j)+(M.Wscaled(nu,k)*mag(k))*base;
    end
else
    G=zeros(4,numel(active));
    for h=1:numel(active)
        [nu,k]=ind2sub(size(M.Wscaled),active(h));map(nu,k)=h;sgn(nu,k)=1;
        base=zeros(4,1);base(nu)=1;
        if strcmp(side,'right'),base=localMul(base,power(:,k));else,base=localMul(power(:,k),base);end
        G(:,h)=(M.Wscaled(nu,k)*mag(k))*base;
    end
end
if any(~isfinite([G(:);value]))
    E.message='Value or zonotope construction overflowed.';return
end
if M.coefficientScale==0
    E.radius=0;E.lower=0;E.upper=0;E.valid=true;E.status='computed';
    E.message='The zero nominal polynomial already vanishes; not a stability assertion.';return
end
Z=qstab.internal.zonotopegauge(G,-value,tol,maxFacets);E.gauge=Z;E.generatorCount=sum(any(G~=0,1));
if Z.infeasible
    missing=all(G==0,2);
    if any(missing & value~=0) && max(abs(value(missing)))<=256*eps*(n+1)*max(amplitude,realmin)
        E.message='Unattainability is within evaluation roundoff; no infinite gauge is asserted.';return
    end
    E.radius=Inf;E.lower=Inf;E.upper=Inf;E.valid=true;E.status='infeasible';
    E.message='The permitted generators cannot cancel the value at this actual quaternion point.';return
end
if ~isfinite(Z.upper)
    E.lower=localScale(Z.lower,M.logScaleRatio);
    E.message='No checked finite primal zonotope representation was found.';return
end
E.radius=localScale(Z.upper,M.logScaleRatio);E.lower=localScale(Z.lower,M.logScaleRatio);
if ~isfinite(E.radius)||(E.radius==0&&Z.upper>0)
    E.radius=NaN;E.lower=0;E.message='Positive component radius is not representable.';return
end
% A primal representation is always reconstructed for validity; Witness only
% controls whether these ordinary data are retained in the public result.
D=zeros(size(M.C));Ccrit=M.C;
for h=1:numel(active)
    [nu,k]=ind2sub(size(M.Wscaled),active(h));j=map(nu,k);
    if j>0,D(nu,k)=(M.coefficientScale*M.Wscaled(nu,k))*sgn(nu,k)*Z.coefficients(j);end
    Ccrit(nu,k)=M.C(nu,k)+D(nu,k);
end
if a==0
    Ccrit=M.C;D(:)=0;
    for nu=1:4
        if M.weights(nu,1)>0,Ccrit(nu,1)=0;D(nu,1)=-M.C(nu,1);end
    end
end
if any(~isfinite(Ccrit(:)))||any(~isfinite(D(:)))
    E.message='Critical polynomial overflowed.';return
end
actual=Ccrit-M.C;ratios=zeros(size(actual));mask=M.weights>0;
ratios(mask)=abs(actual(mask))./M.weights(mask);achieved=max(ratios(:));
fixedOK=all(actual(~mask)==0);res=zeros(4,1);den=0;
for k=1:n+1
    ck=Ccrit(:,k)/M.coefficientScale;
    if strcmp(side,'right'),term=localMul(ck,power(:,k));else,term=localMul(power(:,k),ck);end
    res=res+mag(k)*term;den=den+mag(k)*norm(ck);
end
relativeResidual=norm(res)/max([den,amplitude,realmin]);
valid=fixedOK&&isfinite(achieved)&&relativeResidual<=max(4096*eps*(n+1),10*tol);
if ~valid,E.message='Returned component-box witness failed reconstruction checks.';return;end
E.upper=max(E.radius,achieved)*(1+512*eps*(n+1));
if E.lower>E.upper*(1+100*tol)
    E.lower=0;E.message='Primal and dual values conflict numerically.';return
end
E.lower=min(E.lower,E.upper);E.valid=true;E.status='bounded';
if E.upper-E.lower<=max(4096*eps*(n+1),100*tol)*max(E.upper,realmin),E.status='computed';end
E.message='Pointwise box gauge; this is not minimization over a quaternion similarity class.';
if makeWitness
    E.witness=struct('kind','finiteZero','valid',true,'evaluationSide',side, ...
        'boundaryComponents',q,'classRepresentative',complex(q(1),norm(q(2:4))), ...
        'coefficientComponents',Ccrit,'perturbationComponents',actual, ...
        'polynomial',localParts(Ccrit),'weightedNorm',achieved, ...
        'evaluationResidual',relativeResidual,'fixedComponentsPreserved',fixedOK, ...
        'weights',M.weights,'grade',n,'validated',false);
end
end
function P=localParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function z=localMul(x,y)
z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function y=localScale(x,lr)
if x==0,y=0;elseif isinf(x),y=Inf;else,y=exp(log(x)+lr);end
end
