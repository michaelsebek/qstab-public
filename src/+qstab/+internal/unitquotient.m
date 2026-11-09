function E = unitquotient(M,theta,side,tolerance)
%UNITQUOTIENT Nearest one-quaternion parameter on a complete unit class.
% Left parameter uses RIGHT evaluation; right parameter uses LEFT evaluation.
% The formula is the sphere quotient from qB01, now at cos(theta)+sin(theta)*u.
if ~(isnumeric(theta)&&isreal(theta)&&isscalar(theta)&&isfinite(theta)&&theta>=0&&theta<=pi)
    error('qstab:unitquotient:Angle','Angle must lie in the closed interval [0,pi].');
end
n=M.grade;ks=0:n;
if theta==0,powers=ones(1,n+1);z=1;
elseif theta==pi,powers=(-1).^ks;z=-1;
else,powers=exp(1i*ks*theta);z=complex(cos(theta),sin(theta));end
V0=M.S0*powers.';V1=M.S1*powers.';e0=real(V0);o0=imag(V0);e1=real(V1);o1=imag(V1);
E=struct('valid',false,'radius',NaN,'lower',0,'upper',Inf,'scaledRadius',NaN, ...
    'scaledLower',0,'scaledMinNominal',0,'angle',theta,'boundaryArgument',z, ...
    'parameterComponents',NaN(4,1),'boundaryComponents',NaN(4,1), ...
    'forwardBoundaryComponents',NaN(4,1),'evaluationResidual',NaN, ...
    'witness',[],'message','');
a=e0'*e0+o0'*o0;
if strcmp(side,'left'),cc=mul(conjq(o0),e0);dd=mul(conjq(o1),e1);evaluationSide='right';
else,cc=mul(e0,conjq(o0));dd=mul(e1,conjq(o1));evaluationSide='left';end
c=cc(2:4);d=dd(2:4);cn=norm(c);
amp0=sum(sqrt(sum(M.S0.^2,1)));amp1=sum(sqrt(sum(M.S1.^2,1)));
err0=256*eps*(n+1)*max(amp0,realmin);
m0=min(svd(localMatrix(V0)));E.scaledMinNominal=max(0,m0-err0);
if ~(isfinite(a)&&m0>max(err0,tolerance*sqrt(a)))
    E.message='Nominal separation on this class is unresolved.';return
end
if ~M.hasPerturbation || all([e1;o1]==0)
    E.valid=true;E.radius=Inf;E.lower=Inf;E.upper=Inf;E.scaledRadius=Inf;E.scaledLower=Inf;
    E.message='The perturbation vanishes on this complete boundary class.';return
end
b=e1'*e1+o1'*o1;
% m0^2 avoids cancellation in a-2*norm(c); both formulas are equivalent.
delta0=m0*m0*(a+2*cn);delta01=a*b-4*(c'*d);
if cn>0
    projection=(c'*d)/cn;perp=d-(projection/cn)*c;
    rad=4*((a*projection-b*cn)^2+delta0*(perp'*perp));
else,rad=4*a*a*(d'*d);end
den=delta01+sqrt(max(0,rad));
if ~(isfinite(den)&&den>0&&isfinite(delta0)&&delta0>0)
    E.message='The quotient is numerically singular.';return
end
r2=delta0/den;rs=sqrt(r2);
vec=-(c-r2*d);nv=norm(vec);
if nv<=128*eps*(norm(c)+r2*norm(d)+realmin),vec=[1;0;0];else,vec=vec/nv;end
uq=[0;vec];
if strcmp(side,'left'),f0=e0+mul(o0,uq);f1=e1+mul(o1,uq);
else,f0=e0+mul(uq,o0);f1=e1+mul(uq,o1);end
fscale=max(abs(f1));
if ~(fscale>0&&isfinite(fscale)&&isfinite(rs)&&rs>0)
    E.message='The critical parameter cannot be reconstructed reliably.';return
end
f1s=f1/fscale;invf1=conjq(f1s)/(f1s'*f1s);
if strcmp(side,'left'),ts=-mul(f0,invf1)/fscale;else,ts=-mul(invf1,f0)/fscale;end
t=scalevector(ts,M.logScaleRatio);r=norm(t);
if any(~isfinite(t))||r==0||~isfinite(r)
    E.message='The positive critical parameter exceeds floating-point range.';return
end
% Check the polynomial actually returned in the ORIGINAL variable.
D=zeros(size(M.inputC1));Q=M.inputC0;
for j=1:n+1
    if strcmp(side,'left'),D(:,j)=mul(t,M.inputC1(:,j));else,D(:,j)=mul(M.inputC1(:,j),t);end
    Q(:,j)=Q(:,j)+D(:,j);
end
lambda=[real(z);imag(z)*vec];lambdaInput=lambda;
if strcmp(M.variable,'zinv'),lambdaInput=conjq(lambda);end
if any(~isfinite(Q(:)))||any(~isfinite(D(:)))
    E.message='Critical polynomial coefficients overflow.';return
end
qscale=max(abs(Q(:)));residual=0;value=zeros(4,1);denom=0;pow=[1;0;0;0];
if qscale>0
    for j=1:n+1
        qj=Q(:,j)/qscale;
        if strcmp(side,'left'),value=value+mul(qj,pow);else,value=value+mul(pow,qj);end
        denom=denom+norm(qj);pow=mul(pow,lambdaInput);
    end
    residual=norm(value)/max(denom,realmin);
end
disagreement=abs(norm(ts)-rs)/max(rs,realmin);
if residual>max(1024*eps*(n+1),100*tolerance)||disagreement>100*tolerance
    E.message='The reconstructed polynomial or quotient failed its residual check.';return
end
arithmetic=512*eps*(n+1)*max(1,(amp0+rs*amp1)/m0);
scaledLower=max(0,min(norm(ts),rs)*(1-100*tolerance-arithmetic));
E.valid=true;E.radius=r;E.lower=scalevalue(scaledLower,M.logScaleRatio);
E.upper=r*(1+512*eps*(n+1));E.scaledRadius=rs;E.scaledLower=scaledLower;
E.parameterComponents=t;E.boundaryComponents=lambdaInput;E.forwardBoundaryComponents=lambda;
E.evaluationResidual=residual;
E.witness=struct('kind','finiteZero','valid',true,'side',side,'evaluationSide',evaluationSide, ...
    'variable',M.variable,'testVariable','z','grade',n,'parameterComponents',t, ...
    'parameterNorm',r,'boundaryComponents',lambdaInput,'forwardBoundaryComponents',lambda, ...
    'classRepresentative',complex(real(z),abs(imag(z))), ...
    'coefficientComponents',Q,'perturbationComponents',D,'polynomial',parts(Q), ...
    'evaluationResidual',residual,'formulaWitnessRelativeDifference',disagreement,'validated',false);
E.message='Nearest numerical one-quaternion parameter at this full unit-boundary class.';
end
function Q=localMatrix(v)
Q=[v(1)+1i*v(2),v(3)+1i*v(4);-v(3)+1i*v(4),v(1)-1i*v(2)];
end
function z=mul(x,y)
z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function z=conjq(x),z=[x(1);-x(2:4)];end
function y=scalevalue(x,l)
if x==0,y=0;elseif isinf(x),y=Inf;else,y=exp(log(x)+l);end
end
function y=scalevector(x,l)
y=zeros(size(x));for j=1:numel(x),y(j)=sign(x(j))*scalevalue(abs(x(j)),l);end
end
function P=parts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
