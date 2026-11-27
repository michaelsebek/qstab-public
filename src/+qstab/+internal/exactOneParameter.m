function varargout = exactOneParameter(action,varargin)
%EXACTONEPARAMETER Rational algebra for optional one-parameter certificates.
% All polynomial variables are local symbolic objects. No assumptions,
% preferences, caller variables, paths or random states are changed.
switch action
    case 'require'
        if isempty(which('sym'))
            error('qstab:exact:SymbolicRequired', ...
                ['This optional certificate requires Symbolic Math Toolbox. ', ...
                 'The numerical qstab functions do not require it.']);
        end
    case 'data'
        [varargout{1:nargout}]=readData(varargin{:});
    case 'rational'
        varargout{1}=rational(varargin{:});
    case 'nominal'
        [varargout{1:nargout}]=nominal(varargin{:});
    case 'frequency'
        [varargout{1:nargout}]=frequency(varargin{:});
    case 'ray'
        [varargout{1:nargout}]=ray(varargin{:});
    case 'sturm'
        [varargout{1:nargout}]=sturm(varargin{:});
    case 'count'
        [varargout{1:nargout}]=count(varargin{:});
    case 'isolate'
        [varargout{1:nargout}]=isolate(varargin{:});
    case 'sign'
        varargout{1}=sgn(varargin{1});
    case 'iszero'
        varargout{1}=allZero(varargin{1});
    case 'strings'
        varargout{1}=stringsOf(varargin{1});
    case 'enclosure'
        varargout{1}=doubleEnclosure(varargin{:});
    case 'mul'
        varargout{1}=mul(varargin{:});
    case 'gcd'
        varargout{1}=polyGcd(varargin{:});
    otherwise
        error('qstab:exact:Internal','Unknown exact helper operation.');
end
end

function [A,B,d,source]=readData(P0,P1,d0,R)
A=readPolynomialComponents(P0);B=readPolynomialComponents(P1);
if isempty(d0),d=[];else,d=direction(d0);end
source='exact-binary-input';
if ~isempty(R)
    if ~isstruct(R)||~isscalar(R)||~all(isfield(R,{'P0','P1'}))
        error('qstab:exact:RationalData','RationalData needs P0 and P1 component matrices.');
    end
    unknown=setdiff(fieldnames(R),{'P0','P1','Direction'});
    if ~isempty(unknown),error('qstab:exact:RationalData','Unknown RationalData field %s.',unknown{1});end
    A=override(A,R.P0,'P0');B=override(B,R.P1,'P1');
    if isfield(R,'Direction')
        if isempty(d0),error('qstab:exact:RationalData','Direction is only meaningful for a ray certificate.');end
        if numel(R.Direction)~=4,error('qstab:exact:RationalData','Direction needs four real components.');end
        d=override(d,reshape(R.Direction,4,1),'Direction');
    end
    source='explicit-rational-data';
end
A=trim(A);B=trim(B);
end

function C=readPolynomialComponents(P)
if isa(P,'sym')
    if isempty(P)||~isvector(P),error('qstab:exact:Coefficients','A symbolic coefficient input must be a nonempty vector.');end
    P=reshape(P,1,numel(P)); C=[real(P);imag(P);sym(zeros(2,numel(P)))];
    C=rational(C);
elseif isstruct(P)&&isscalar(P)&&all(isfield(P,{'a','b','c','d'}))&& ...
        (isa(P.a,'sym')||isa(P.b,'sym')||isa(P.c,'sym')||isa(P.d,'sym'))
    names={'a','b','c','d'};C=[];n=[];
    for k=1:4
        t=P.(names{k});
        if isempty(t)||~isvector(t)||(~isempty(n)&&numel(t)~=n)
            error('qstab:exact:Coefficients','Symbolic component arrays must be equal nonempty vectors.');
        end
        n=numel(t);C=[C;reshape(rational(t),1,n)]; %#ok<AGROW>
    end
else
    C=rational(qstab.polycomponents(P));
end
end

function d=direction(u)
if isa(u,'sym')
    if ~isscalar(u),error('qstab:exact:Direction','The quaternion direction must be scalar.');end
    d=rational([real(u);imag(u);sym(0);sym(0)]);
elseif isstruct(u)&&isscalar(u)&&all(isfield(u,{'a','b','c','d'}))&& ...
        any(structfun(@(v)isa(v,'sym'),u))
    d=rational([u.a;u.b;u.c;u.d]);
    if numel(d)~=4,error('qstab:exact:Direction','The quaternion direction must be scalar.');end
else
    [a,b,c,e]=qstab.components(u);
    if numel(a)~=1,error('qstab:exact:Direction','The quaternion direction must be scalar.');end
    d=rational([a;b;c;e]);
end
if allZero(d),error('qstab:exact:Direction','The quaternion direction must be nonzero.');end
end

function S=override(old,new,name)
if ~ismatrix(new)||size(new,1)~=4
    error('qstab:exact:RationalData','RationalData.%s must be a two-dimensional 4-by-N component matrix.',name);
end
S=rational(new);
if ~strcmp(name,'Direction')&&ismatrix(S)&&size(S,1)==4
    S=trim(S);old=trim(old);
end
if ~isequal(size(S),size(old))||size(S,1)~=4
        error('qstab:exact:RationalData','RationalData.%s must have the original 4-by-N component shape (apart from exact trailing zeros).',name);
end
% An explicit rational source may round to the supplied doubles, but may
% not silently certify another family or discard a nonzero component.
if ~isequal(double(S),double(old))||~isequal(zeroMask(S),zeroMask(old))
    error('qstab:exact:RationalMismatch','RationalData.%s does not round to the supplied coefficients.',name);
end
end

function S=rational(x)
if isa(x,'sym')
    S=x;
    requireRationals(S);return
elseif isnumeric(x)||islogical(x)
    if ~isreal(x)||any(~isfinite(x(:)))
        error('qstab:exact:Rational','Exact component data must be finite and real.');
    end
    S=sym(double(x),'f');return
elseif ischar(x)||(isstring(x)&&isscalar(x))
    if ischar(x)&&~isrow(x)
        error('qstab:exact:Rational','Each rational text value must be one row of characters.');
    end
    s=strtrim(char(x));
    atom='[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?';
    if isempty(regexp(s,['^',atom,'(?:\s*/\s*',atom,')?$'],'once'))
        error('qstab:exact:Rational','Use an integer, decimal or fraction for exact rational text.');
    end
    slash=strfind(s,'/');
    if isempty(slash)
        S=decimalAtom(s);
    else
        numerator=decimalAtom(strtrim(s(1:slash-1)));
        denominator=decimalAtom(strtrim(s(slash+1:end)));
        if isequal(denominator,sym(0))
            error('qstab:exact:Rational','An exact rational denominator must be nonzero.');
        end
        S=numerator/denominator;
    end
    return
elseif iscell(x)||isstring(x)
    S=sym(zeros(size(x)));
    for k=1:numel(x)
        if iscell(x),v=x{k};else,v=x(k);end
        if ~(isscalar(v)||ischar(v)),error('qstab:exact:Rational','Each rational component must be scalar.');end
        S(k)=rational(v);
    end
    return
else
    error('qstab:exact:Rational','Use numeric, symbolic rational or rational-text data.');
end
end

function q=decimalAtom(s)
% Parse the decimal spelling as integers and powers of ten. Passing decimal
% text directly to SYM can create symbolic floating-point atoms in MATLAB.
% Such atoms must never enter a polynomial GCD or a Sturm certificate.
negative=false;
if s(1)=='-'||s(1)=='+'
    negative=s(1)=='-';s=s(2:end);
end
exponent=sym(0);at=regexp(s,'[eE]','once');
if ~isempty(at)
    exponent=integerText(s(at+1:end));s=s(1:at-1);
end
dot=strfind(s,'.');
if ~isempty(dot)
    exponent=exponent-sym(numel(s)-dot);s(dot)=[];
end
q=integerText(s);
if ~isequal(q,sym(0)),q=q*sym(10)^exponent;end
if negative,q=-q;end
end

function q=integerText(s)
% Only a sign and decimal integer digits reach this helper.
negative=false;
if s(1)=='-'||s(1)=='+'
    negative=s(1)=='-';s=s(2:end);
end
s=regexprep(s,'^0+','');
if isempty(s),q=sym(0);else,q=sym(s);end
if negative,q=-q;end
end

function requireRationals(S)
% A symbolic approximation is not exact input, even when its numerical
% value happens to be integral. Do not rationalize VPA data silently.
if exist('OCTAVE_VERSION','builtin')==0
    accepted=isSymType(S,'integer | rational');
    ok=all(accepted(:));
else
    % Conservative compatibility for auxiliary Octave runs, which do not
    % implement ISSYMTYPE. Only canonical integer numerator/denominator
    % spellings are accepted; floating-point spellings are rejected.
    [n,d]=numden(S);ok=true;
    for k=1:numel(S)
        ok=ok&&~isempty(regexp(char(n(k)),'^-?[0-9]+$','once'))&& ...
            ~isempty(regexp(char(d(k)),'^-?[0-9]+$','once'))&& ...
            ~isequal(d(k),sym(0));
        if ~ok,break;end
    end
end
if ~ok
    error('qstab:exact:Rational', ...
        ['Symbolic components must be exact integers or rational numbers. ', ...
         'VPA approximations, variables and nonrational constants are not exact input; ', ...
         'use decimal or fraction text for exact decimal source data.']);
end
end

function C=trim(C)
while size(C,2)>1&&allZero(C(:,end)),C(:,end)=[];end
end

function [ok,P]=nominal(C)
s=sym('qstabExactS');n=size(C,2)-1;
z=sym(zeros(4,1));
for j=1:size(C,2),z=z+C(:,j)*s^(j-1);end
N=expand(sum(z.^2));a=fliplr(coeffs(N,s,'All'));
D=sym(zeros(1,2*n));H=hurwitz(a,2*n,2*n);
ok=n>=1;
for j=1:2*n
    D(j)=expand(det(H(1:j,1:j)));
    ok=ok&&sgn(D(j))>0;
end
P=struct('strictlyHurwitz',logical(ok),'method','exact-Hurwitz-determinants-of-real-norm', ...
    'degree',n,'normCoefficients',{stringsOf(a)},'hurwitzDeterminants',{stringsOf(D)});
end

function H=hurwitz(a,degree,m)
H=sym(zeros(m));
for i=1:m
    for j=1:m
        k=degree-2*j+i;
        if k>=0&&k<numel(a),H(i,j)=a(k+1);end
    end
end
end

function [n0,n1,delta0,delta1,delta01,x]=frequency(A,B,side)
x=sym('qstabExactX');
[n0,c0]=normParts(A,x,side);[n1,c1]=normParts(B,x,side);
delta0=expand(n0^2-4*x*sum(c0.^2));
delta1=expand(n1^2-4*x*sum(c1.^2));
delta01=expand(n0*n1-4*x*sum(c0.*c1));
end

function [n,c]=normParts(A,x,side)
e=sym(zeros(4,1));o=e;
for k=0:size(A,2)-1
    if mod(k,2)==0,e=e+A(:,k+1)*(-x)^(k/2);
    else,o=o+A(:,k+1)*(-x)^((k-1)/2);end
end
n=expand(sum(e.^2)+x*sum(o.^2));
oc=[o(1);-o(2:4)];
if strcmp(side,'left'),c=mul(oc,e);else,c=mul(e,oc);end
c=expand(c(2:4));
end

function [h,N,g,s,C]=ray(A,B,d,side)
g=sym('qstabExactGamma');s=sym('qstabExactS');
n=size(A,2)-1;B(:,end+1:n+1)=sym(0);C=A;
for k=1:n+1
    if strcmp(side,'left'),v=mul(d,B(:,k));else,v=mul(B(:,k),d);end
    C(:,k)=A(:,k)+g*v;
end
p=sym(zeros(4,1));for k=0:n,p=p+C(:,k+1)*s^k;end
N=expand(sum(p.^2));a=fliplr(coeffs(N,s,'All'));
h=expand(simplify(det(hurwitz(a,2*n,2*n-1))));
end

function z=mul(a,b)
z=[a(1)*b(1)-sum(a(2:4).*b(2:4)); ...
    a(1)*b(2)+b(1)*a(2)+a(3)*b(4)-a(4)*b(3); ...
    a(1)*b(3)+b(1)*a(3)+a(4)*b(2)-a(2)*b(4); ...
    a(1)*b(4)+b(1)*a(4)+a(2)*b(3)-a(3)*b(2)];
end

function [S,f]=sturm(p,x)
p=expand(p);
% Refuse approximate coefficients before any exact GCD/root-count claim.
requireRationals(coeffs(p,x,'All'));
if allZero(p),error('qstab:exact:Internal','Sturm input must not vanish identically.');end
d=diff(p,x);
if allZero(d),S={p};f=p;return;end
c=polyGcd(p,d,x);[f,rem]=quorem(p,c,x);
if ~allZero(expand(rem)),error('qstab:exact:Internal','Square-free division failed.');end
f=expand(f);S={normalize(f,x),normalize(diff(f,x),x)};
while true
    [~,r]=quorem(S{end-1},S{end},x);r=expand(-r);
    if allZero(r),break;end
    S{end+1}=normalize(r,x); %#ok<AGROW>
end
end

function g=polyGcd(a,b,x)
a=expand(a);b=expand(b);
while ~allZero(b)
    [~,r]=quorem(a,b,x);a=b;b=expand(r);
end
if allZero(a),g=sym(0);else
    c=coeffs(a,x,'All');g=expand(a/c(1));
end
end

function p=normalize(p,x)
c=coeffs(p,x,'All');scale=abs(c(1));p=expand(p/scale);
end

function [n,vlo,vhi]=count(S,x,lo,hi)
vlo=variations(S,x,lo);vhi=variations(S,x,hi);n=vlo-vhi;
end

function v=variations(S,x,a)
signs=zeros(1,numel(S));
for k=1:numel(S)
    if isnumeric(a)&&isinf(a)
        c=coeffs(S{k},x,'All');signs(k)=sgn(c(1));
        if a<0,signs(k)=signs(k)*(-1)^(numel(c)-1);end
    else
        signs(k)=sgn(subs(S{k},x,a));
    end
end
signs(signs==0)=[];v=sum(signs(1:end-1).*signs(2:end)<0);
end

function [lo,hi,I]=isolate(p,x,tol,maxSteps)
[S,f]=sturm(p,x);zero=sym(0);[total,v0,vi]=count(S,x,zero,Inf);
I=struct('positiveRootCount',total,'isolated',false,'isRationalRoot',false, ...
    'iterations',0,'widthTargetMet',false,'variationAtZero',v0,'variationAtInfinity',vi, ...
    'squareFreePolynomial',char(f),'squareFreeCoefficients',{stringsOf(fliplr(coeffs(f,x,'All')))}, ...
    'sturmSequence',{{}},'rootsInEnclosure',0);
I.sturmSequence=cellfun(@char,S,'UniformOutput',false);
if total==0,lo=Inf;hi=Inf;I.isolated=true;I.widthTargetMet=true;return;end
c=coeffs(f,x,'All');bound=sym(1);
for k=2:numel(c)
    a=abs(c(k)/c(1));if sgn(a+1-bound)>0,bound=a+1;end
end
lo=zero;hi=sym(1);while sgn(hi-bound)<=0,hi=2*hi;end
target=rational(tol);vlo=v0;vhi=vi;
for k=1:maxSteps
    mid=(lo+hi)/2;vmid=variations(S,x,mid);nr=v0-vmid;
    if nr==1&&sgn(subs(f,x,mid))==0
        lo=mid;hi=mid;I.isolated=true;I.isRationalRoot=true;I.iterations=k;break
    end
    if nr>0,hi=mid;vhi=vmid;else,lo=mid;vlo=vmid;end
    nr=vlo-vhi;I.iterations=k;
    if nr==1&&sgn(hi-lo-target)<=0,I.isolated=true;break;end
end
if sgn(hi-lo)==0,nr=1;else,nr=vlo-vhi;end
I.rootsInEnclosure=nr;I.isolated=I.isolated||nr==1;
I.widthTargetMet=sgn(hi-lo-target)<=0;
end

function y=sgn(x)
x=simplify(x);
if ~isscalar(x)||~isempty(symvar(x))
    error('qstab:exact:Undecidable','An exact scalar sign requires a constant, not a free polynomial variable.');
end
if allZero(x),y=0;
elseif isAlways(x>0,'Unknown','false'),y=1;
elseif isAlways(x<0,'Unknown','false'),y=-1;
else,error('qstab:exact:Undecidable','An exact rational sign could not be decided.');end
end

function tf=allZero(x)
mask=zeroMask(x);tf=all(mask(:));
end

function mask=zeroMask(x)
% This is a polynomial/rational-function identity test, not a question
% whether a parameter-dependent expression happens to vanish at a point.
% Exact expansion gives a canonical zero numerator without ISALWAYS's
% TruthUnknown warnings on nonzero polynomials.
[n,~]=numden(x);n=expand(n);mask=false(size(x));
for k=1:numel(x),mask(k)=isequal(n(k),sym(0));end
end

function C=stringsOf(S)
C=cell(size(S));for k=1:numel(S),C{k}=char(S(k));end
end

function b=doubleEnclosure(lo,hi)
if isnumeric(lo)&&isinf(lo),b=[lo hi];return;end
b=[double(lo),double(hi)];
if isinf(b(1))&&b(1)>0,b(1)=realmax;end
if isfinite(b(1))&&sgn(rational(b(1))-lo)>0,b(1)=b(1)-eps(b(1));end
if isfinite(b(2))&&sgn(rational(b(2))-hi)<0,b(2)=b(2)+eps(b(2));end
end
