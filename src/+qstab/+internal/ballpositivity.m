function S = ballpositivity(M,r,maxCells,maxDepth,maxTime)
%BALLPOSITIVITY Sufficient CT lower-radius proof on the whole compact interval.
% For eta=1,2,Inf, H(t)=B(t)'*B(t)-r^2*rho(t)^2*I is polynomial,
% or a finite list of such polynomials for eta=1 (max of monomials).
% Bernstein MATRIX controls are subdivided until all are positive definite.
% This gives a lower bound, never a counterexample when the check fails.
% Uses ordinary double precision and explicit guards, not outward rounding.
S=struct('attempted',false,'accepted',false,'radius',r,'processedCells',0, ...
    'acceptedCells',zeros(0,4),'matrixCount',0,'validated',false, ...
    'reason','The exponent is not one of 1,2,Inf.');
if ~any(M.exponent==[1 2 Inf]),return;end
S.attempted=true;
if r==0,S.accepted=true;S.reason='Nonnegativity gives the zero lower bound.';return;end
if maxCells<1||maxTime<=0,S.reason='No shortcut budget.';return;end
n=M.grade;d=2*n;t0=tic;
rwork=exp(log(r)-M.logScaleRatio);
if ~isfinite(rwork)||rwork==0,S.reason='Target scaling is not representable.';return;end
c1=ones(1,n+1);c2=ones(1,d+1);
for k=1:n,c1(k+1)=c1(k)*(n-k+1)/k;end
for k=1:d,c2(k+1)=c2(k)*(d-k+1)/k;end
B=zeros(2,2,n+1);
for k=0:n
    v=M.Cscaled(:,k+1)*1i^k/c1(k+1);
    B(:,:,k+1)=[v(1)+1i*v(2),v(3)+1i*v(4);-v(3)+1i*v(4),v(1)-1i*v(2)];
end
HH=zeros(2,2,d+1);aa=zeros(1,d+1);gpoly=zeros(1,d+1);
for i=0:n
    for j=0:n
        factor=c1(i+1)*c1(j+1)/c2(i+j+1);
        HH(:,:,i+j+1)=HH(:,:,i+j+1)+factor*(B(:,:,i+1)'*B(:,:,j+1));
        aa(i+j+1)=aa(i+j+1)+factor*norm(B(:,:,i+1),2)*norm(B(:,:,j+1),2);
        if isinf(M.exponent)
            gpoly(i+j+1)=gpoly(i+j+1)+M.weightsScaled(i+1)*M.weightsScaled(j+1)/c2(i+j+1);
        end
    end
end
if M.exponent==2
    for k=0:n,gpoly(2*k+1)=M.weightsScaled(k+1)^2/c2(2*k+1);end
end
if M.exponent==1,indices=find(M.weightsScaled>0);else,indices=0;end
S.matrixCount=numel(indices);allH=zeros(4*numel(indices),d+1);scale=0;
for z=1:numel(indices)
    if M.exponent==1
        gp=zeros(1,d+1);k=indices(z)-1;gp(2*k+1)=M.weightsScaled(k+1)^2/c2(2*k+1);
    else,gp=gpoly;end
    for k=1:d+1
        H=HH(:,:,k)-rwork^2*gp(k)*eye(2);
        H=(H+H')/2;allH(4*z-3:4*z,k)=H(:);
        scale=max(scale,aa(k)+rwork^2*gp(k));
    end
end
if any(~isfinite(allH(:)))||~isfinite(scale)
    S.reason='Bernstein product arithmetic overflowed.';return
end
baseGuard=1024*eps*(n+1)^2*max(scale,realmin);
queue=struct('left',0,'right',1,'depth',0,'controls',allH);
S.reason='Shortcut budget or positivity test did not resolve the full domain.';
while ~isempty(queue)
    if S.processedCells>=maxCells||toc(t0)>=maxTime,return;end
    node=queue(end);queue(end)=[];S.processedCells=S.processedCells+1;
    smallest=Inf;
    for k=1:d+1
        for j=1:numel(indices)
            A=reshape(node.controls(4*j-3:4*j,k),2,2);
            smallest=min(smallest,min(real(builtin('eig',(A+A')/2))));
        end
    end
    guard=baseGuard*(node.depth+1);
    if smallest>guard
        S.acceptedCells(end+1,:)=[node.left,node.right,node.depth,smallest-guard]; %#ok<AGROW>
        continue
    end
    if node.depth>=maxDepth,return;end
    middle=(node.left+node.right)/2;
    if middle==node.left||middle==node.right,return;end
    [L,R]=split(node.controls);
    queue(end+1)=struct('left',middle,'right',node.right,'depth',node.depth+1,'controls',R); %#ok<AGROW>
    queue(end+1)=struct('left',node.left,'right',middle,'depth',node.depth+1,'controls',L); %#ok<AGROW>
end
S.accepted=true;S.reason='Guarded positive Bernstein matrix controls cover the entire compact CT boundary.';
end
function [L,R]=split(C)
n=size(C,2);L=zeros(size(C));R=L;D=C;L(:,1)=D(:,1);R(:,n)=D(:,n);
for k=2:n,D=(D(:,1:end-1)+D(:,2:end))/2;L(:,k)=D(:,1);R(:,n-k+1)=D(:,end);end
end
