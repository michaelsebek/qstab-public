function R = zonotopegauge(G,x,tolerance,maxFacets)
%ZONOTOPEGAUGE Gauge of G*[-1,1]^N in its intrinsic REAL dimension (<=4).
% Dual support directions give lower bounds. A checked primal representation
% gives an upper bound. Finite facet enumeration can be budgeted; incomplete
% enumeration or uncertain rank NEVER becomes a fabricated infinity.
% This is floating point, not outward-rounded interval verification.
x=double(x(:));G=double(G);m=size(G,1);N=size(G,2);
R=struct('status','undetermined','radius',NaN,'lower',0,'upper',Inf, ...
    'coefficients',zeros(N,1),'numericalRank',0,'spanResidual',NaN, ...
    'primalResidual',Inf,'facetCount',0,'complete',false,'dualDirection',zeros(m,1), ...
    'valid',false,'infeasible',false,'validated',false);
if isempty(x)||m~=numel(x)||m>4||~isreal(G)||~isreal(x)||any(~isfinite([G(:);x]))
    error('qstab:box:GaugeInput','Gauge inputs must be finite real arrays with one to four rows.');
end
if all(x==0)
    R.status='computed';R.radius=0;R.upper=0;R.primalResidual=0;R.spanResidual=0;
    R.complete=true;R.valid=true;return
end
if N==0||all(G(:)==0)
    R=localInfinite(R);return
end
% An exactly zero coordinate row cannot change that coordinate of x.
if any(all(G==0,2) & x~=0)
    R=localInfinite(R);return
end
gs=max(abs(G(:)));xs=max(abs(x));A=G/gs;b=x/xs;
S=localSolve(A,b,tolerance,maxFacets,0);
R.numericalRank=S.rank;R.spanResidual=S.spanResidual;R.primalResidual=S.residual;
R.facetCount=S.facets;R.complete=S.complete;R.dualDirection=S.dual;
ratio=log(xs)-log(gs);
R.lower=localScale(S.lower,ratio);R.upper=localScale(S.upper,ratio);
R.coefficients=S.coefficients*localScale(1,ratio);
if any(~isfinite(R.coefficients))||isnan(R.lower)||isnan(R.upper)
    R.coefficients=zeros(N,1);R.upper=Inf;R.status='undetermined';return
end
if isfinite(R.upper)
    R.radius=R.upper;
    if R.lower>R.upper*(1+100*tolerance)
        R.lower=0;R.status='undetermined';return
    end
    R.lower=min(R.lower,R.upper);
    R.valid=true;R.status='bounded';
    if R.upper-R.lower<=max(2048*eps,100*tolerance)*max(R.upper,realmin)
        R.status='computed';
    end
end
end
function R=localInfinite(R)
R.status='infeasible';R.radius=Inf;R.lower=Inf;R.upper=Inf;
R.infeasible=true;R.valid=true;R.complete=true;
end
function S=localSolve(A,b,tol,budget,depth)
[m,N]=size(A);S=struct('lower',0,'upper',Inf,'coefficients',zeros(N,1), ...
    'rank',0,'spanResidual',Inf,'residual',Inf,'facets',0,'complete',false,'dual',zeros(m,1));
if all(b==0),S.upper=0;S.spanResidual=0;S.residual=0;S.complete=true;return;end
if N==0||all(A(:)==0),return;end
active=find(any(A~=0,1));B=A(:,active);
[U,D,V]=svd(B,'econ');s=diag(D);s=s(1:min(size(B)));
threshold=64*eps*max(size(B))*max(s);d=sum(s>threshold);S.rank=d;
if d==0,return;end
Ur=U(:,1:d);br=Ur'*b;span=b-Ur*br;
S.spanResidual=norm(span)/max(norm(b),realmin);
% Every direction is a valid support inequality, regardless of numerical rank.
for j=1:size(U,2)
    c=U(:,j);[bound,~]=localSupport(B,b,c);if bound>S.lower,S.lower=bound;S.dual=c;end
end
if norm(span)>0
    c=span/norm(span);[bound,~]=localSupport(B,b,c);
    if bound>S.lower,S.lower=bound;S.dual=c;end
end
% A minimum-two-norm feasible solution supplies a conservative primal upper.
coef=V(:,1:d)*(br./s(1:d));[ok,res,up]=localCheck(B,b,coef,tol);
if ok,S.upper=up;S.coefficients(active)=coef;S.residual=res;end
if depth>=4||budget==0,return;end
H=Ur'*B;k=d-1;nn=size(B,2);best=-Inf;bestC=[];
if k==0,idx=zeros(1,0);else,idx=1:k;end
more=true;
while more && S.facets<budget
    if k==0
        c=Ur(:,1);independent=true;
    else
        [~,T,Z]=svd(H(:,idx)');ts=diag(T);ts=ts(1:k);
        independent=min(ts)>128*eps*max(1,norm(H(:,idx),'fro'));
        if independent,c=Ur*Z(:,end);else,c=[];end
    end
    if independent
        S.facets=S.facets+1;
        [bound,raw]=localSupport(B,b,c);
        if bound>S.lower,S.lower=bound;S.dual=c;end
        if raw>best,best=raw;bestC=c;end
    else
        % Count rejected tuples too: the budget is a work bound.
        S.facets=S.facets+1;
    end
    [idx,more]=localNext(idx,nn);
end
S.complete=~more;
if isempty(bestC)||~isfinite(best)||best<0,return;end
c=bestC;if c'*b<0,c=-c;end
dots=c'*B;support=sum(abs(dots));
if support==0,return;end
tau=(c'*b)/support;
% Exposed-face reconstruction. Nearly orthogonal columns are left free and
% subsequently checked in the ORIGINAL equality, never discarded as zero.
free=abs(dots)<=256*eps*max(1,norm(B,'fro'));
f=zeros(nn,1);f(~free)=tau*sign(dots(~free));residual=b-B*f;
if any(free)
    sub=localSolve(B(:,free),residual,tol,max(0,budget-S.facets),depth+1);
    S.facets=S.facets+sub.facets;
    if isfinite(sub.upper),f(free)=sub.coefficients;end
end
[ok,res,up]=localCheck(B,b,f,tol);
if ok&&up<S.upper
    S.upper=up;S.coefficients(:)=0;S.coefficients(active)=f;S.residual=res;
end
end
function [v,raw]=localSupport(A,b,c)
numerator=abs(c'*b);denominator=sum(abs(c'*A));
errN=128*eps*sum(abs(c).*abs(b));errD=128*eps*sum(abs(c)'*abs(A));
v=max(0,numerator-errN)/max(denominator+errD,realmin);
if denominator==0,raw=-Inf;else,raw=numerator/denominator;end
end
function [ok,res,up]=localCheck(A,b,c,tol)
res=norm(A*c-b)/max(norm(b),realmin);
up=max(abs(c));ok=all(isfinite(c))&&isfinite(res)&&res<=max(4096*eps,10*tol);
if ~ok,up=Inf;end
end
function [idx,more]=localNext(idx,N)
k=numel(idx);more=false;if k==0,return;end
for j=k:-1:1
    if idx(j)<N-k+j
        idx(j)=idx(j)+1;
        for h=j+1:k,idx(h)=idx(h-1)+1;end
        more=true;return
    end
end
end
function y=localScale(x,lr)
if x==0,y=0;elseif isinf(x),y=Inf;else,y=exp(log(x)+lr);end
end
