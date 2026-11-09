function M = qstab_eigen_match(actual,reference,tolerance)
%QSTAB_EIGEN_MATCH Bottleneck bijection of canonical complex multisets.
% Multiplicities are never removed. Row j in permutation is the actual index
% assigned to reference j. Distances are pairwise scale-relative, with the
% largest reference modulus as a floor (one only for an all-zero reference).
% Invalid actual output is diagnosed, not silently truncated or projected.
if nargin<3,tolerance=1e-7;end
if ~isnumeric(tolerance)||~isreal(tolerance)||~isscalar(tolerance)|| ...
        ~isfinite(tolerance)||tolerance<=0
    error('qstab:eigencomparison:Tolerance','Tolerance must be finite and positive.');
end
if ~isnumeric(reference)||(~isempty(reference)&&~isvector(reference))||any(~isfinite(reference(:)))
    error('qstab:eigencomparison:Reference','The reference must be a finite numeric vector.');
end
y=complex(real(reference(:)),abs(imag(reference(:))));
M=struct('status','INVALID','ok',false,'error',Inf,'permutation',[], ...
    'reference',y,'actual',[],'tolerance',tolerance,'scale',1,'message','');
if ~isnumeric(actual)||(~isempty(actual)&&~isvector(actual))||any(~isfinite(actual(:)))
    M.message='Actual representatives must be a finite numeric vector.';return
end
x=complex(real(actual(:)),abs(imag(actual(:))));M.actual=x;
if numel(x)~=numel(y)
    M.status='MISMATCH';M.message='Algebraic counts differ.';return
end
n=numel(x);
if n==0,M.status='MATCHED';M.ok=true;M.error=0;M.message='Both multisets are empty.';return;end
base=max(abs(y));if base==0,base=1;end;M.scale=base;
s=max(base,bsxfun(@max,abs(x),abs(y.')));
% Divide BEFORE subtraction, avoiding overflow under common large scales.
C=abs(bsxfun(@rdivide,x,s)-bsxfun(@rdivide,y.',s));
levels=unique(C(:));low=1;high=numel(levels);
while low<high
    middle=floor((low+high)/2);
    if perfect(C<=levels(middle)),high=middle;else,low=middle+1;end
end
M.error=levels(low);[~,M.permutation]=perfect(C<=M.error);
M.ok=M.error<=tolerance;
if M.ok,M.status='MATCHED';M.message='A multiplicity-preserving bijection meets tolerance.';
else,M.status='MISMATCH';M.message='No multiplicity-preserving bijection meets tolerance.';end
end
function [yes,match]=perfect(E)
n=size(E,1);match=zeros(1,n);yes=true;
for row=1:n
    [found,match,~]=augment(row,E,match,false(1,n));
    if ~found,yes=false;return;end
end
end
function [ok,match,seen]=augment(row,E,match,seen)
ok=false;
for col=find(E(row,:))
    if seen(col),continue;end
    seen(col)=true;
    if match(col)==0,match(col)=row;ok=true;return;end
    [found,proposal,seen]=augment(match(col),E,match,seen);
    if found,match=proposal;match(col)=row;ok=true;return;end
end
end
