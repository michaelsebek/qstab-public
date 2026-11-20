function [gradient,Hessian,Info] = rootmotion(family,lambda,varargin)
%ROOTMOTION Local motion of a simple algebraic zero of an affine family.
%   [G,H,Info] = qstab.rootmotion(F,LAMBDA) differentiates the selected
%   simple zero of the REAL NORM polynomial of a scalar quaternion family.
%   F is a qstab.family definition of kind affine, interval, or bounded
%   oneparameter. LAMBDA is its finite complex algebraic representative.
%   The default expansion point is the midpoint of F's real parameter box.
%
%   [G,H,Info] = qstab.rootmotion(F,LAMBDA,'Parameters',THETA) uses the
%   supplied real parameter column THETA. Parameter coordinates and order
%   are exactly those returned by qstab.familycomponents(F). In particular,
%   interval coordinates are coefficient COMPONENTS, not normalized widths.
%   THETA must belong to the declared closed parameter box.
%
%   G is the complex first derivative column and H the complex symmetric
%   Hessian. For a real parameter displacement h,
%       lambda(theta+h) = lambda + G.'*h + h.'*H*h/2 + O(norm(h)^3).
%   At a stationary imaginary-axis contact the quadratic real-part form
%   used in the interval article is Q=-real(H)/2. This is local sensitivity,
%   not a certificate for a finite parameter box or a multiple zero.
%
%   Name/value options:
%     Parameters []       Midpoint, or a finite real vector of parameters.
%     Tolerance  1e-10    Relative root-residual and derivative guards.
%   Info.valid is false and G,H are NaN when the root is not verified as
%   simple. Info retains residuals, the derivative separation, and the full
%   parameter convention. No eigensolver or finite differences are used.
%   Norm products use ordinary squares, not absolute squares.
%
%   Example (standard MATLAB data, coefficients ascending):
%     F = qstab.family('affine',[1 2 1],{[0 1]},[-.1 .1]);
%     [g,H,I] = qstab.rootmotion(F,-1); % repeated norm zero: unresolved
%
%   See also qstab.family, qstab.normpoly, qstab.roots, qstab.isstable.
try
ip=inputParser;ip.FunctionName='qstab.rootmotion';ip.PartialMatching=false;
addParameter(ip,'Parameters',[]);
addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<1);
parse(ip,varargin{:});o=ip.Results;
if ~isnumeric(lambda)||~isscalar(lambda)||~isfinite(lambda)
    error('qstab:rootmotion:Root','LAMBDA must be a finite complex numeric scalar.');
end
F=qstab.familycomponents(family);
if strcmp(F.kind,'polytope')
    error('qstab:rootmotion:Family','Use an affine family with explicit independent real coordinates.');
end
m=size(F.bounds,1);
theta=o.Parameters;
if isempty(theta),theta=.5*F.bounds(:,1)+.5*F.bounds(:,2);end
if ~isnumeric(theta)||~isreal(theta)||numel(theta)~=m||any(~isfinite(theta(:)))
    error('qstab:rootmotion:Parameters','Parameters must be a finite real vector with one entry per affine coordinate.');
end
theta=double(theta(:));
if any(theta<F.bounds(:,1)|theta>F.bounds(:,2))
    error('qstab:rootmotion:Parameters','The expansion point must belong to the declared parameter box.');
end
C=F.base;
for k=1:m,C=C+theta(k)*F.generators(:,:,k);end
if any(~isfinite(C(:)))
    error('qstab:rootmotion:Range','The polynomial at the expansion point is not representable.');
end
scale=max(abs([C(:);F.generators(:)]));
gradient=complex(nan(m,1));Hessian=complex(nan(m,m));
Info=struct('schema','qstab.rootmotion.1','status','undetermined','valid',false, ...
    'lambda',lambda,'parameters',theta,'parameterNames',{F.parameterNames}, ...
    'relativeRootResidual',Inf,'derivativeSeparation',0,'normValue',NaN, ...
    'normDerivative',NaN,'normSecondDerivative',NaN,'coefficientScale',scale, ...
    'message','','validated',false,'scope','localSimpleNormZero', ...
    'assurance','floating-point implicit differentiation; no finite-box certificate');
if scale==0,Info.message='The zero norm polynomial has no simple zero.';return;end
C=C/scale;G=F.generators./scale;
[p,dp,ddp]=jet(C,lambda);
V=complex(zeros(4,m));DV=V;
for k=1:m,[V(:,k),DV(:,k)]=jet(G(:,:,k),lambda);end
N=p.'*p;Nz=2*(p.'*dp);Nzz=2*(dp.'*dp+p.'*ddp);
den=sum(abs(p).^2);sep=2*sum(abs(p).*abs(dp));
Info.normValue=N;Info.normDerivative=Nz;Info.normSecondDerivative=Nzz;
Info.relativeRootResidual=abs(N)/max(den,realmin);
Info.derivativeSeparation=abs(Nz)/max(sep,realmin);
if any(~isfinite([p;dp;ddp;V(:);DV(:);N;Nz;Nzz]))
    Info.message='Polynomial jets overflowed.';return;
end
if Info.relativeRootResidual>o.Tolerance
    Info.message='The supplied representative does not pass the norm-root residual guard.';return;
end
if Nz==0||Info.derivativeSeparation<=o.Tolerance
    Info.message='A simple norm zero is not numerically separated from a multiple zero.';return;
end
Nt=2*(V.'*p);Nzt=2*(V.'*dp+DV.'*p);Ntt=2*(V.'*V);
gradient=-Nt/Nz;
Hessian=-(Ntt+Nzt*gradient.'+gradient*Nzt.'+Nzz*(gradient*gradient.'))/Nz;
Info.valid=all(isfinite([gradient(:);Hessian(:)]));
if Info.valid,Info.status='computed';Info.message='Local simple-root derivatives computed analytically.';
else,Info.message='Root derivatives are not representable.';gradient(:)=NaN;Hessian(:)=NaN;end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'),throwAsCaller(problem);end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
function [p,dp,ddp]=jet(C,z)
p=C(:,end);dp=zeros(4,1);ddp=zeros(4,1);
for k=size(C,2)-1:-1:1
    ddp=ddp*z+2*dp;dp=dp*z+p;p=p*z+C(:,k);
end
end
