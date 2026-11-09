function [C,S] = sampleFamily(F,N,o)
%SAMPLEFAMILY Draw whole admissible polynomial members before mapping them.
%   C(:,r+1,n) is the [real;i;j;k] coefficient of v^r in sample n, in
%   the INPUT variable (no delay reversal). Each row of S.parameters or
%   S.weights identifies the complete member C(:,:,n).
method=o.sampling; seed=o.seed; window=o.box;
if ~((ischar(method)&&isrow(method)) || (isstring(method)&&isscalar(method)&&~ismissing(method)))
    error('qstab:sample:Sampling','Sampling must be grid or random.');
end
method=lower(strtrim(char(method)));
if ~any(strcmp(method,{'grid','random'}))
    error('qstab:sample:Sampling','Sampling must be grid or random.');
end
if ~isequal(seed,[]) && ~(isnumeric(seed)&&isreal(seed)&&isscalar(seed)&& ...
        isfinite(seed)&&seed>=0&&seed<=2^32-1&&seed==fix(seed))
    error('qstab:sample:Seed','Seed must be [] or an integer from 0 to 2^32-1.');
end
[D,v,g,period,unrestricted]=qstab.internal.familyPlotModel(F,o);
if ~isequal(window,[]) && ~strcmp(D.kind,'oneparameter')
    error('qstab:sample:SampleBounds','A sampling box applies only to a oneparameter family.');
end
bounds=D.bounds; originalBounds=bounds;
if unrestricted, originalBounds=zeros(0,2); end
if strcmp(D.kind,'oneparameter')
    if isequal(window,[]) && isempty(bounds)
        error('qstab:sample:SampleBounds', ...
            'An unrestricted quaternion has no finite uniform grid or uniform distribution. Supply a finite 4-by-2 sampling box.');
    elseif ~isequal(window,[])
        if ~(isnumeric(window)&&isreal(window)&&isequal(size(window),[4,2])&& ...
                all(isfinite(window(:)))&&all(window(:,1)<=window(:,2)))
            error('qstab:sample:SampleBounds','Use four finite rows [lower upper] in real,i,j,k order.');
        end
        window=full(double(window));
        if ~isempty(bounds) && any(window(:,1)<bounds(:,1) | window(:,2)>bounds(:,2))
            error('qstab:sample:SampleBounds','The sampling box must lie inside the specified family bounds.');
        end
        bounds=window;
    end
end
% Cap the largest requested data array before allocating plotted samples.
if strcmp(D.kind,'polytope'), wd=size(D.vertices,3); else, wd=size(bounds,1); end
if N*max([1,4*(g+1),wd])>1e7
    error('qstab:sample:Count','Requested sample exceeds 1e7 coefficient or parameter entries. Reduce N or Grade.');
end
draw=@rand;
if strcmp(method,'random') && ~isempty(seed)
    if exist('RandStream','class')==8
        stream=RandStream('mt19937ar','Seed',double(seed));
        draw=@(varargin)rand(stream,varargin{:});
    else
        old=rng; restore=onCleanup(@()rng(old)); %#ok<NASGU>
        rng(double(seed),'twister');
    end
end
T=zeros(N,0); W=zeros(N,0); order=[];
if strcmp(D.kind,'polytope')
    m=size(D.vertices,3);
    if m==1
        W=ones(N,1); design='singleton';
    elseif strcmp(method,'random')
        E=-log(max(realmin,draw(N,m)));       % independent Exp(1) coordinates
        W=bsxfun(@rdivide,E,sum(E,2));       % uniform on the weight simplex
        design='uniform simplex weights';
    else
        [W,design,order]=simplexGrid(N,m);
    end
    C=reshape(reshape(D.vertices,4*(D.degree+1),m)*W.',4,D.degree+1,N);
else
    free=bounds(:,1)<bounds(:,2);
    m=nnz(free); T=repmat(bounds(:,1).',N,1);
    if m==0
        U=zeros(N,0); design='singleton';
    elseif strcmp(method,'random')
        U=draw(N,m); design='uniform independent free coordinates';
    else
        [U,design,order]=boxGrid(N,m);
    end
    if m>0, T(:,free)=scaleBox(U,bounds(free,:)); end
    if strcmp(D.kind,'interval')
        % For a component box its parameters ARE the varying coefficient
        % entries. Avoid N products with every sparse unit-direction page.
        flat=repmat(D.base(:),1,N);
        if ~isempty(T)
            idx=sub2ind([4,D.degree+1],D.parameterIndex(:,1),1+D.parameterIndex(:,2));
            flat(idx,:)=T.';
        end
        C=reshape(flat,4,D.degree+1,N);
    else
        C=repmat(D.base,1,1,N);
        for j=1:size(T,2)
            C=C+bsxfun(@times,D.generators(:,:,j),reshape(T(:,j),1,1,N));
        end
    end
end
if g>D.degree, C(:,end+1:g+1,:)=0; end
if any(~isfinite(C(:)))
    error('qstab:sample:CoefficientRange','A sampled polynomial is not representable in finite double coefficients.');
end
S=struct('kind',D.kind,'grade',g,'variable',v,'samplingPeriod',period, ...
    'count',N,'sampling',method,'seed',seed,'design',design,'gridOrder',order, ...
    'parameters',T,'weights',W,'bounds',bounds,'originalBounds',originalBounds,'side',D.side, ...
    'parameterNames',{D.parameterNames},'parameterIndex',zeros(0,2),'isCertificate',false);
if isfield(D,'parameterIndex'),S.parameterIndex=D.parameterIndex;end
end

function T=scaleBox(U,B)
% Avoid overflow in hi-lo for opposite signs; keep fixed coordinates exact.
T=zeros(size(U));
for j=1:size(U,2)
    a=B(j,1); b=B(j,2); t=U(:,j);
    if a<0 && b>0
        col=(1-t)*a+t*b;
    else
        col=a+t*(b-a);
    end
    col(t==0)=a; col(t==1)=b;
    T(:,j)=min(b,max(a,col)); % protect against a 1-ulp endpoint overshoot
end
end

function [U,design,q]=boxGrid(N,m)
% q = levels in the largest complete equal tensor grid fitting N.
if N==1, U=.5*ones(1,m); design='box midpoint'; q=1; return; end
if m==1, U=linspace(0,1,N).'; design='uniform interval grid'; q=N; return; end
q=floor(N^(1/m));
while (q+1)^m<=N, q=q+1; end
while q^m>N, q=q-1; end
U=zeros(N,m);
if q>=2
    used=q^m; n=(0:used-1).';
    for j=1:m
        U(1:used,j)=mod(n,q)/(q-1);
        n=floor(n/q);
    end
    design='complete Cartesian grid'; h=q-1;
else
    U(2,:)=1; used=2; h=1;
    design='opposite corners'; % a complete corner mesh would exceed N
end
if used<N
    if mod(h,2)==1
        used=used+1; U(used,:)=.5;
        design=[design ' + midpoint'];
    end
    if used<N
        U(used+1:end,:)=interiorGrid(N-used,m,2*h);
        design=[design ' + stratified grid'];
    end
end
end

function [W,design,h]=simplexGrid(N,m)
% h = denominator of the largest complete barycentric lattice fitting N.
% Enumerate ONLY accepted compositions; never construct an Ngrid > budget.
if N==1, W=ones(1,m)/m; design='simplex centroid'; h=0; return; end
if m==2
    w=linspace(0,1,N).'; W=[1-w w]; design='uniform segment grid';h=N-1;return
end
if N<m
    W=zeros(N,m);
    if N==2, ix=1; else, ix=round(linspace(1,m,N-1)); end
    W(sub2ind([N,m],1:N-1,ix))=1; W(N,:)=1/m;
    design='selected simplex vertices + centroid';h=0;return
end
lo=1;hi=N;
while hi>lo
    mid=lo+ceil((hi-lo)/2);
    if latticeCount(mid,m,N)<=N,lo=mid;else,hi=mid-1;end
end
h=lo; used=latticeCount(h,m,N);
W=zeros(N,m); a=zeros(1,m); a(1)=h;
for n=1:used
    W(n,:)=a/h;
    if n<used
        j=find(a(1:end-1)>0,1,'last');
        tail=a(end);           % all between j and the last position are zero
        a(j)=a(j)-1; a(end)=0; a(j+1)=tail+1;
    end
end
design='complete barycentric grid';
if used<N && mod(h,m)~=0
    used=used+1; W(used,:)=1/m; design=[design ' + centroid'];
end
if used<N
    r=N-used;
    cuts=sort(interiorGrid(r,m-1,h*m),2);
    W(used+1:end,:)=diff([zeros(r,1),cuts,ones(r,1)],1,2);
    design=[design ' + stratified simplex grid'];
end
end

function n=latticeCount(h,m,limit)
% (h+m-1 choose m-1), stopped as soon as the result exceeds the budget.
k=min(h,m-1); n=1;
for j=1:k
    n=round(n*((h+m-1-k+j)/j));
    if n>limit,n=limit+1;return;end
end
end

function X=interiorGrid(N,m,avoid)
% Deterministic Latin grid: one point in each of N strata per coordinate.
% The first offset cannot hit a multiple of 1/avoid. Distinct offsets for
% other axes keep every barycentric spacing positive after sorting cuts.
X=zeros(N,m); k=(0:N-1).'; offset=0.5/avoid;
for j=1:m
    if j==1 || N==1, perm=k;
    else
        step=1+mod(floor(N*mod(j*.6180339887498949,1)),N);
        while gcd(step,N)~=1, step=mod(step,N)+1; end
        shift=floor(N*mod(j*.4142135623730950,1));
        perm=mod(k*step+shift,N);
    end
    X(:,j)=(perm+(offset+(j-1)/(2*m)))/N;
end
end
