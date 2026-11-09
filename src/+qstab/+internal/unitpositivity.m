function S = unitpositivity(M,r,side,maxCells,maxDepth,maxTime)
%UNITPOSITIVITY Sufficient whole-semicircle one-parameter radius lower bound.
% z=(1-x+i*x)/(1-x-i*x), x in [0,1]. The common degree-n denominator
% cancels in the ratio. Left parameter: B0'*B0-r^2*B1'*B1 > 0;
% right parameter: B0*B0'-r^2*B1*B1' > 0. Bernstein MATRIX controls
% retain the same quaternion direction in numerator and denominator.
S=struct('attempted',false,'accepted',false,'radius',r,'processedCells',0, ...
    'acceptedCells',zeros(0,4),'validated',false,'method','unit-circle-matrix-Bernstein', ...
    'parameterization','z=(1-x+i*x)/(1-x-i*x), x in [0,1]', ...
    'reason','No shortcut budget.');
if maxCells<1||maxTime<=0,return;end
S.attempted=true;
if r==0,S.accepted=true;S.reason='Nonnegativity supplies the zero lower bound.';return;end
rw=exp(log(r)-M.logScaleRatio);
if ~isfinite(rw)||rw==0,S.reason='Target scaling is unrepresentable.';return;end
n=M.grade;d=2*n;timer=tic;common=max(1,rw);
C0=M.S0/common;C1=M.S1*(rw/common);
if any(M.S0(:)~=0&C0(:)==0)||any(M.S1(:)~=0&C1(:)==0)
    S.reason='Target scaling would erase coefficient data.';return
end
[B0,amp0]=controls(C0);[B1,amp1]=controls(C1);
bn=binomials(n);bd=binomials(d);H=zeros(2,2,d+1);absolute=zeros(1,d+1);
for i=0:n
    for j=0:n
        f=bn(i+1)*bn(j+1)/bd(i+j+1);
        if strcmp(side,'left')
            term=B0(:,:,i+1)'*B0(:,:,j+1)-B1(:,:,i+1)'*B1(:,:,j+1);
        else
            term=B0(:,:,i+1)*B0(:,:,j+1)'-B1(:,:,i+1)*B1(:,:,j+1)';
        end
        H(:,:,i+j+1)=H(:,:,i+j+1)+f*term;
        absolute(i+j+1)=absolute(i+j+1)+f*(amp0(i+1)*amp0(j+1)+amp1(i+1)*amp1(j+1));
    end
end
for k=1:d+1,H(:,:,k)=(H(:,:,k)+H(:,:,k)')/2;end
baseGuard=2048*eps*(n+1)^3*max([absolute,realmin]);
if any(~isfinite(H(:)))||~isfinite(baseGuard)
    S.reason='Matrix-control arithmetic is unrepresentable.';return
end
queue=struct('left',0,'right',1,'depth',0,'controls',reshape(H,4,d+1));
S.reason='The sufficient positivity test did not finish; use boundary subdivision.';
while ~isempty(queue)
    if S.processedCells>=maxCells||toc(timer)>=maxTime,return;end
    node=queue(end);queue(end)=[];S.processedCells=S.processedCells+1;
    smallest=Inf;
    for k=1:d+1
        A=reshape(node.controls(:,k),2,2);
        smallest=min(smallest,min(real(builtin('eig',(A+A')/2))));
    end
    if smallest>baseGuard*(1+node.depth)
        S.acceptedCells(end+1,:)=[node.left,node.right,node.depth,smallest-baseGuard*(1+node.depth)]; %#ok<AGROW>
        continue
    end
    if node.depth>=maxDepth,return;end
    splitPoint=node.left+(node.right-node.left)/2;
    if splitPoint==node.left||splitPoint==node.right,return;end
    [leftControls,rightControls]=split(node.controls);
    queue(end+1)=struct('left',splitPoint,'right',node.right,'depth',node.depth+1,'controls',rightControls); %#ok<AGROW>
    queue(end+1)=struct('left',node.left,'right',splitPoint,'depth',node.depth+1,'controls',leftControls); %#ok<AGROW>
end
S.accepted=true;S.reason='Positive guarded Bernstein matrix controls cover the full upper unit semicircle.';
end
function [B,amp]=controls(C)
n=size(C,2)-1;B=zeros(2,2,n+1);amp=zeros(1,n+1);
for k=0:n
    b=1;
    for j=1:k,b=multiplyLinear(b,[1,1i]);end
    for j=1:n-k,b=multiplyLinear(b,[1,-1i]);end
    c=C(:,k+1);Q=[c(1)+1i*c(2),c(3)+1i*c(4);-c(3)+1i*c(4),c(1)-1i*c(2)];
    for j=1:n+1
        B(:,:,j)=B(:,:,j)+b(j)*Q;
        % Absolute Bernstein product terms are bounded by one before summing.
        amp(j)=amp(j)+norm(c);
    end
end
end
function c=multiplyLinear(b,a)
m=numel(b)-1;c=zeros(1,m+2);
for j=0:m
    c(j+1)=c(j+1)+(m+1-j)/(m+1)*b(j+1)*a(1);
    c(j+2)=c(j+2)+(j+1)/(m+1)*b(j+1)*a(2);
end
end
function b=binomials(n)
b=ones(1,n+1);for j=1:n,b(j+1)=b(j)*(n-j+1)/j;end
end
function [L,R]=split(C)
n=size(C,2);L=zeros(size(C));R=L;D=C;L(:,1)=D(:,1);R(:,n)=D(:,n);
for k=2:n,D=.5*D(:,1:end-1)+.5*D(:,2:end);L(:,k)=D(:,1);R(:,n-k+1)=D(:,end);end
end
