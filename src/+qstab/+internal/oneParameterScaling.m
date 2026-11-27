function [C0,C1,S] = oneParameterScaling(C0,C1,requested)
%ONEPARAMETERSCALING Reversible real rescaling of a continuous polynomial law.
%   s=tau*sigma and Pj(s)=2^e*PjWorking(sigma). The SAME positive
%   coefficient factor is applied to both polynomials: the quaternion gain,
%   its multiplication side and every ray parameter remain unchanged.
if nargin<3,requested='auto';end
automatic=(ischar(requested)||(isstring(requested)&&isscalar(requested))) ...
    &&strcmpi(char(requested),'auto');
if ~automatic && ~(isnumeric(requested)&&isreal(requested)&&isscalar(requested) ...
        &&isfinite(requested)&&requested>0)
    error('qstab:oneparameter:FrequencyScale','FrequencyScale must be ''auto'' or a positive finite real scalar.');
end
n=size(C0,2)-1;
tau=1;
if automatic
    peak=max(abs(C0(:)));lead=max(abs(C0(:,end)));constant=max(abs(C0(:,1)));
    if n>0&&peak>0&&lead>0&&constant>0 && ...
            (lead/peak<sqrt(eps)||constant/peak<sqrt(eps))
        exponent=round((log2(constant)-log2(lead))/n);
        exponent=max(-1022,min(1023,exponent));
        tau=pow2(1,exponent);
    end
else
    tau=double(requested);
end
% A power-of-two common scale avoids overflow in coefficient powers. For
% ordinary inputs and tau=1 retain the original coefficients bit for bit.
allC=[C0,C1];peak=max(abs(allC(:)));exponent=0;
if tau~=1 || peak>pow2(1,256) || (peak>0&&peak<pow2(1,-256))
    exponents=[weightedlogs(C0,tau),weightedlogs(C1,tau)];
    finite=exponents(isfinite(exponents));
    if ~isempty(finite),exponent=floor(max(finite));end
end
S=struct('schema','qstab.oneparameter-scaling.1','automatic',automatic, ...
    'frequencyScale',tau,'coefficientScalePowerOfTwo',exponent, ...
    'applied',tau~=1||exponent~=0,'originalVariable','s','workingVariable','sigma', ...
    'mapping','s=frequencyScale*sigma; Pj(s)=2^coefficientScalePowerOfTwo*PjWorking(sigma); t unchanged', ...
    'parameterUnchanged',true,'degreePreserved',true,'assurance','floating-point');
if S.applied
    try
        D0=transform(C0,tau,exponent);D1=transform(C1,tau,exponent);
    catch problem
        if tau==1&&strcmp(problem.identifier,'qstab:oneparameter:ScalingRange')
            % An optional common amplitude scale must not erase an input
            % when the parameter scale itself lies outside double range.
            % The existing margin range guards can report that case.
            S.coefficientScalePowerOfTwo=0;S.applied=false;return
        end
        rethrow(problem)
    end
    C0=D0;C1=D1;
end
end

function v=weightedlogs(C,tau)
v=log2(abs(C));
v=bsxfun(@plus,v,(0:size(C,2)-1)*log2(tau));
v=v(:).';
end

function D=transform(C,tau,e)
[mt,et]=log2(tau);[m,b]=log2(C);k=0:size(C,2)-1;
% mt lies in [1/2,1). Splitting its power also prevents underflow for
% unusually high grades, without discarding a nonzero coefficient.
z=k*log2(mt);zi=floor(z);zf=2.^(z-zi);
if all(k<512),zf=mt.^k;zi=zeros(size(k));end
D=pow2(bsxfun(@times,m,zf),bsxfun(@plus,b,et*k+zi-e));
if any(~isfinite(D(:)))||any((C(:)~=0)&(D(:)==0))
    error('qstab:oneparameter:ScalingRange', ...
        'The requested frequency scaling cannot preserve every nonzero coefficient in double precision. Choose a different FrequencyScale.');
end
end
