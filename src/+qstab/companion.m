function [A,B,info] = companion(P,varargin)
%COMPANION Complex companion matrix or generalized companion pencil.
%   [A,B,INFO] = qstab.companion(P) constructs a pencil A-lambda*B for
%   P(v)=P0+P1*v+...+Pg*v^g, using the complex adjoint of each coefficient.
%   P may be a scalar coefficient vector/tensor, a square matrix coefficient
%   tensor, a vector cell array {P0,...,Pg}, or a component struct a,b,c,d.
%   Coefficients are ascending and act from the left. The pencil preserves
%   infinite eigenvalues when the declared leading coefficient is singular.
%   A = qstab.companion(P) instead returns the normalized companion matrix;
%   this one-output form requires a numerically nonsingular leading block.
%   ...,'Grade',G pads with zero high coefficients without silently trimming.
%   ...,'Variable','s'|'z'|'zinv' sets metadata only: no reciprocal transform
%   is performed here. All spectral points are in the actual input variable.
%
%   Examples (ascending coefficients):
%   [A,B,I] = qstab.companion([2 3 1]);
%   e = eig(A,B);  % full complex embedded spectrum, not n quaternion classes
%   disp(e);
%
%   Interpretation:
%   A and B are COMPLEX embedded matrices of order twice the quaternion companion order. They are not quaternion companion matrices and must not be wrapped as such. One output normalizes the leading block; two outputs preserve the generalized pencil.
%
%   Reference: docs/reference/companion.html.
%   See also qstab.about.

grade = []; variable = 's';
if mod(numel(varargin),2)~=0
    error('qstab:InvalidOptions','Options must be name/value pairs.');
end
for k=1:2:numel(varargin)
    name=lower(char(varargin{k}));
    switch name
        case 'grade', grade=varargin{k+1};
        case 'variable', variable=lower(char(varargin{k+1}));
        otherwise, error('qstab:UnknownOption','Unknown option: %s.',name);
    end
end
if ~ismember(variable,{'s','z','zinv'})
    error('qstab:InvalidVariable','Variable must be s, z or zinv.');
end
[C,shape,inputGrade,actualDegree] = localCoefficients(P);
if isempty(grade), grade=inputGrade; end
if ~(isnumeric(grade) && isreal(grade) && isscalar(grade) && ...
        ~isnan(grade) && ~isinf(grade) && grade>=inputGrade && grade==fix(grade))
    error('qstab:InvalidGrade','Grade must be an integer no smaller than the supplied grade.');
end
n=shape(1); m=2*n;
for k=inputGrade+2:grade+1
    C{k}=zeros(m);
end
cn=zeros(1,grade+1);
for k=1:grade+1, cn(k)=norm(C{k},'fro'); end
if n==0
    error('qstab:EmptyPolynomial','Polynomial coefficients must be nonempty.');
end
info=struct('grade',double(grade),'inputGrade',inputGrade, ...
    'actualDegree',actualDegree,'quaternionSize',shape, ...
    'quaternionOrder',n,'isScalar',n==1,'variable',variable, ...
    'coefficientOrder','ascending','coefficientNorms',cn, ...
    'leadingRcond',rcond(C{end}),'constantRcond',rcond(C{1}), ...
    'leadingIsZero',all(C{end}(:)==0),'constantIsZero',all(C{1}(:)==0), ...
    'isZero',all(cn==0),'embeddingResidual',0,'assurance','floating-point');
if nargout<2 && info.leadingRcond<=100*eps
    error('qstab:SingularLeadingCoefficient', ...
        'Use [A,B,info]=qstab.companion(P) for a singular or ill-conditioned leading coefficient.');
end
if grade==0
    A=zeros(0); B=zeros(0); return
end
A=zeros(m*grade); B=eye(m*grade);
for j=1:grade
    A(1:m,(j-1)*m+(1:m))=-C{grade-j+1};
end
if grade>1
    A(m+1:end,1:m*(grade-1))=eye(m*(grade-1));
end
B(1:m,1:m)=C{end};
if nargout<2
    for j=1:grade
        A(1:m,(j-1)*m+(1:m))=-(C{end}\C{grade-j+1});
    end
    if any(isnan(A(:)) | isinf(A(:)))
        error('qstab:ArithmeticOverflow','Normalized companion overflowed; use the pencil form.');
    end
end
end

function [C,shape,g,degree] = localCoefficients(P)
if iscell(P)
    if isempty(P) || ~isvector(P)
        error('qstab:InvalidPolynomial','Coefficient cells must form a nonempty vector.');
    end
    g=numel(P)-1; C=cell(1,g+1); shape=[]; degree=-Inf;
    for k=1:g+1
        [X,xi]=qstab.chi(P{k});
        if isempty(shape), shape=xi.quaternionSize; end
        if ~isequal(shape,xi.quaternionSize)
            error('qstab:CoefficientSizeMismatch','All coefficient matrices must have the same size.');
        end
        C{k}=X;
        if any(X(:)~=0), degree=k-1; end
    end
else
    [a,b,c,d]=qstab.components(P);
    if isempty(a) || ndims(a)>3
        error('qstab:InvalidPolynomial','Use a nonempty coefficient vector or a three-dimensional coefficient tensor.');
    end
    if isvector(a) || isscalar(a)
        a=reshape(a,1,1,[]); b=reshape(b,1,1,[]);
        c=reshape(c,1,1,[]); d=reshape(d,1,1,[]);
    end
    shape=[size(a,1),size(a,2)]; g=size(a,3)-1;
    C=cell(1,g+1); degree=-Inf;
    for k=1:g+1
        C{k}=qstab.chi(struct('a',a(:,:,k),'b',b(:,:,k), ...
            'c',c(:,:,k),'d',d(:,:,k)));
        if any(C{k}(:)~=0), degree=k-1; end
    end
end
if shape(1)~=shape(2)
    error('qstab:NotSquarePolynomial','The initial companion implementation requires square coefficient matrices.');
end
end
