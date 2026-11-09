function E = ballpoint(M,z,makeWitness,tolerance)
%BALLPOINT Minimum weighted perturbation on one algebraic class.
% Complex z is only a class representative, not a chosen quaternion slice.
z=complex(real(z),abs(imag(z)));n=M.grade;k=0:n;a=abs(z);
E=struct('radius',NaN,'lower',0,'upper',Inf,'valid',false, ...
    'classRepresentative',z,'scaledMinimumModulus',NaN,'scaledEnvelope',NaN, ...
    'powerScaleLog',0,'nearZero',false,'witness',[], 'message','');
if ~isfinite(a)
    E.message='Argument norm is outside the representable range.';return
end
if a==0
    w=zeros(1,n+1);w(1)=1;phase=ones(1,n+1);
elseif a<=1
    w=a.^k;phase=(z/a).^k;
else
    w=(1/a).^(n:-1:0);phase=(z/a).^k;E.powerScaleLog=n*log(a);
end
% Do not reinterpret underflowed nonzero monomials as fixed coefficients.
if a>0&&any(w==0 & (any(M.Cscaled~=0,1)|M.weightsScaled>0))
    E.message='A nonzero monomial underflowed at this argument.';return
end
terms=bsxfun(@times,M.Cscaled,w.*phase);
if any(any(M.Cscaled~=0 & bsxfun(@gt,w,zeros(size(M.Cscaled))) & terms==0))
    E.message='A nonzero coefficient contribution underflowed.';return
end
v=sum(terms,2);e=real(v);o=imag(v);
B=localMatrix(v);s=svd(B);m=s(end);
g=M.weightsScaled.*w;
if any(g==0 & M.weightsScaled>0 & w>0)
    E.message='A nonzero uncertainty sensitivity underflowed.';return
end
rho=qstab.internal.ballnorm(g,M.dualExponent);
cnorms=zeros(1,n+1);for jj=1:n+1,cnorms(jj)=norm(M.Cscaled(:,jj));end
amp=sum(cnorms.*w);
err=128*eps*(n+1)*max(amp,realmin);
rerr=128*eps*(n+1)*rho;
E.scaledMinimumModulus=m;E.scaledEnvelope=rho;E.nearZero=m<=err;
if any(~isfinite([v(:);m;rho;err]))
    E.message='Evaluation overflowed.';return
end
if M.coefficientScale==0
    E.radius=0;E.lower=0;E.upper=0;E.valid=true;
    E.message='The centre is identically zero; this is not a stable nominal family.';
    return
end
if rho==0
    E.valid=true;
    if m>err
        E.radius=Inf;E.lower=Inf;E.upper=Inf;
        E.message='No allowed coefficient changes the value at this class.';
    else
        E.radius=NaN;E.lower=0;E.upper=Inf;
        E.message='Zero sensitivity at a numerically singular centre; membership is unresolved.';
    end
    return
end
E.radius=localScale(m/rho,M.logScaleRatio);
E.lower=localScale(max(0,m-err)/(rho+rerr),M.logScaleRatio);
if rho>rerr
    E.upper=localScale((m+err)/(rho-rerr),M.logScaleRatio);
end
if ~isfinite(E.radius)||(E.radius==0&&m>0)||isnan(E.lower)||isnan(E.upper)
    E.radius=NaN;E.lower=0;E.upper=Inf;
    E.message='Positive radius overflows or underflows double precision.';return
end
E.valid=true;E.message='Numerical pointwise distance to the algebraic zero class; not a frequency minimum.';
if ~makeWitness,return;end
% For right evaluation |e+o*u|^2=D+2 Im(conj(o)*e).Im(u).
c=localMul([o(1);-o(2:4)],e);cv=c(2:4);cn=norm(cv);
if cn==0,u=[1;0;0];else,u=-cv/cn;end
lambda=[real(z);imag(z)*u];
if a==0,unit=[1;0;0;0];else,unit=lambda/a;end
qvalue=e+localMul(o,[0;u]);target=norm(qvalue);
if abs(target-m)>max(512*eps*(n+1),100*tolerance)*max(amp,realmin)
    E.valid=false;E.message='Class minimizer disagrees with the complex singular value.';return
end
y=zeros(1,n+1);active=g>0;
if target>0
    if M.exponent==1
        [gmax,idx]=max(g);y(idx)=target/gmax;
    elseif isinf(M.exponent)
        y(active)=target/rho;
    else
        y(active)=(target/rho)*(g(active)/rho).^(M.dualExponent-1);
    end
    vunit=-qvalue/target;
else
    vunit=zeros(4,1);
end
amplitudes=M.weightsScaled.*y;
delta=zeros(4,n+1);qpower=[1;0;0;0];cu=[unit(1);-unit(2:4)];
for j=1:n+1
    delta(:,j)=amplitudes(j)*localMul(vunit,qpower);
    qpower=localMul(qpower,cu);
end
% Work from the exact supplied centre, not from a renormalized copy.
physicalDelta=delta*M.coefficientScale;Ccrit=M.C;
for j=find(M.weights>0)
    Ccrit(:,j)=M.C(:,j)+physicalDelta(:,j);
end
if a==0
    % At the real zero class only q_0 contributes.  rho>0 above guarantees
    % that this coefficient is free.  Set the EXACT cancellation directly,
    % rather than normalize q_0 and multiply it back (which can leave a
    % roundoff-sized constant whose relative evaluation residual is one).
    % All higher coefficients remain bit-for-bit equal to the input.
    physicalDelta=zeros(size(M.C));physicalDelta(:,1)=-M.C(:,1);
    Ccrit=M.C;Ccrit(:,1)=0;
end
if any(~isfinite(Ccrit(:)))||any(~isfinite(physicalDelta(:)))
    E.valid=false;E.message='Reconstructed coefficients overflowed.';return
end
actualDelta=Ccrit-M.C;
changedNorm=zeros(1,n+1);
for jj=1:n+1,changedNorm(jj)=norm(actualDelta(:,jj)/M.coefficientScale);end
ratios=zeros(1,n+1);mask=M.weightsScaled>0;
ratios(mask)=changedNorm(mask)./M.weightsScaled(mask);
wn=localScale(qstab.internal.ballnorm(ratios,M.exponent),M.logScaleRatio);
% Check the actual returned polynomial, in homogeneous powers to avoid overflow.
qpower=[1;0;0;0];res=zeros(4,1);normalizer=0;
for j=1:n+1
    cj=Ccrit(:,j)/M.coefficientScale;
    res=res+w(j)*localMul(cj,qpower);normalizer=normalizer+w(j)*norm(cj);
    qpower=localMul(qpower,unit);
end
residual=norm(res)/max(normalizer,realmin);
fixedOK=all(all(actualDelta(:,~mask)==0));
witnessOK=isfinite(wn)&&fixedOK&&residual<=max(1024*eps*(n+1),100*tolerance);
E.witness=struct('kind','finiteZero','valid',witnessOK, ...
    'boundaryComponents',lambda,'classRepresentative',z, ...
    'coefficientComponents',Ccrit,'perturbationComponents',actualDelta, ...
    'polynomial',localParts(Ccrit),'weightedNorm',wn, ...
    'evaluationResidual',residual,'fixedCoefficientsPreserved',fixedOK, ...
    'directionAgreement',abs(target-m)/max(amp,realmin), ...
    'weights',M.weights,'exponent',M.exponent,'grade',n,'validated',false);
if ~witnessOK
    E.valid=false;E.message='Reconstructed finite-class witness failed its numerical residual check.';
else
    % Keep the rounding guard as well as the achieved perturbation norm.
    E.upper=max(E.upper,wn*(1+256*eps*(n+1)));
end
end
function B=localMatrix(v)
B=[v(1)+1i*v(2),v(3)+1i*v(4);-v(3)+1i*v(4),v(1)-1i*v(2)];
end
function z=localMul(x,y)
z=[x(1)*y(1)-x(2:4)'*y(2:4);x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function y=localScale(x,logScale)
if x==0,y=0;elseif isinf(x),y=Inf;else,y=exp(log(x)+logScale);end
end
function P=localParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
