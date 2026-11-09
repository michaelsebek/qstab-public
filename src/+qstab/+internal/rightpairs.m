function [Z,lambda,info] = rightpairs(a,b,c,d)
%RIGHTPAIRS Numeric complex-adjoint core for RIGHT quaternion eigenpairs.
% Internal real-component input; Z is the first n columns of chi(V):
% V=U+W*j corresponds to Z=[U;-conj(W)]. No private toolbox is called.
% This body is also embedded, byte-identically, as a local function in the
% separate replacement qeig.m. Neither distribution depends on the other.
% Close computed classes are grouped only for vector SELECTION; their
% separate paired eigenvalues and algebraic counts are not collapsed.

n=size(a,1);
scale=max([abs(a(:));abs(b(:));abs(c(:));abs(d(:));0]);
if scale==0,scale=1;end
as=a/scale;bs=b/scale;cs=c/scale;ds=d/scale;
if any(a(:)~=0 & as(:)==0)||any(b(:)~=0 & bs(:)==0)|| ...
        any(c(:)~=0 & cs(:)==0)||any(d(:)~=0 & ds(:)==0)
    error('rightpairs:ScalingUnderflow', ...
        'Common matrix scaling would erase nonzero components; no spectrum is returned.');
end
X=[complex(as,bs),complex(cs,ds);-complex(cs,-ds),complex(as,-bs)];
xnorm=norm(X,'fro');
pairTolerance=1e-8;
residualTolerance=1e-8;
clusterTolerance=256*eps*max(1,xnorm);
selectionTolerance=32*eps*max(1,2*n);
if n==0
    W=complex(zeros(0));mu=complex(zeros(0,1));
else
    [W,T]=builtin('eig',X);
    mu=diag(T);
end
if any(~isfinite(mu))||any(~isfinite(W(:)))
    error('rightpairs:NonfiniteSolver','The embedded eigensolver returned nonfinite data.');
end
% Pair by FULL complex conjugacy, never by adjacency after real-part sorting.
left=(1:2*n).';pairs=zeros(n,2);lambdaScaled=complex(zeros(n,1));
pairResidual=zeros(n,1);
for k=1:n
    [~,j]=max(abs(imag(mu(left))));first=left(j);left(j)=[];
    [gap,j]=min(abs(mu(left)-conj(mu(first))));second=left(j);left(j)=[];
    z=mu(first)/2+conj(mu(second))/2;
    lambdaScaled(k)=complex(real(z),abs(imag(z)));
    pairs(k,:)=[first second];pairResidual(k)=gap/max(realmin,xnorm);
end
[~,order]=sortrows([real(lambdaScaled),imag(lambdaScaled)],[1 2]);
lambdaScaled=lambdaScaled(order);pairs=pairs(order,:);pairResidual=pairResidual(order);
% Group only indistinguishable computed representatives, at an eps-scale
% threshold. Distinct entries remain distinct in lambdaScaled and D.
remaining=(1:n).';groups=cell(0,1);
while ~isempty(remaining)
    seed=remaining(1);
    take=abs(lambdaScaled(remaining)-lambdaScaled(seed))<=clusterTolerance;
    groups{end+1,1}=remaining(take); %#ok<AGROW>
    remaining(take)=[];
end
Z=complex(zeros(2*n,n));selected=zeros(n,1);reflected=false(n,1);
novelty=zeros(n,1);clusterIndex=zeros(n,1);
for g=1:numel(groups)
    targets=groups{g};ids=reshape(pairs(targets,:).',[],1);
    pool=W(:,ids);flipped=false(numel(ids),1);
    for j=1:numel(ids)
        % chi(A)*J*conj(w)=conj(mu)*J*conj(w); transform the vector TOO.
        if imag(mu(ids(j)))<0
            pool(:,j)=rightpairs_anti(pool(:,j));flipped(j)=true;
        end
        vnorm=norm(pool(:,j));
        if vnorm==0||~isfinite(vnorm)
            error('rightpairs:ZeroVector','The embedded eigensolver returned an invalid column.');
        end
        pool(:,j)=pool(:,j)/vnorm;
    end
    available=(1:numel(ids)).';B=complex(zeros(2*n,0));
    for j=1:numel(targets)
        scores=zeros(numel(available),1);
        for t=1:numel(available)
            r=pool(:,available(t));
            % Reorthogonalization is for SELECTION ONLY. The returned
            % eigenvector is never mixed with other eigenvalues.
            r=r-B*(B'*r);r=r-B*(B'*r);
            scores(t)=norm(r);
        end
        [best,t]=max(scores);pick=available(t);available(t)=[];
        dest=targets(j);Z(:,dest)=pool(:,pick);
        selected(dest)=ids(pick);reflected(dest)=flipped(pick);
        novelty(dest)=best;clusterIndex(dest)=g;
        % span_H{v_1,...} is span_C{z_1,J*conj(z_1),...}.
        % This prevents two complex columns representing the same H-line
        % from displacing an independent quaternionic eigenvector.
        candidates=[B,pool(:,pick),rightpairs_anti(pool(:,pick))];
        [U,S,~]=svd(candidates,'econ');s=diag(S);
        B=U(:,s>selectionTolerance*max(1,s(1)));
    end
end
lambda=complex(real(lambdaScaled)*scale,imag(lambdaScaled)*scale);
outputFinite=all(isfinite(lambda));
outputUnderflow=any(real(lambdaScaled)~=0 & real(lambda)==0)|| ...
    any(imag(lambdaScaled)~=0 & imag(lambda)==0);
returnedScaled=complex(real(lambda)/scale,imag(lambda)/scale);
colResidual=zeros(n,1);
for k=1:n
    colResidual(k)=norm(X*Z(:,k)-Z(:,k)*returnedScaled(k))/ ...
        max(realmin,(xnorm+abs(returnedScaled(k)))*norm(Z(:,k)));
end
if n==0
    globalResidual=0;basisRank=0;rankTolerance=0;rcondBasis=1;
    singular=zeros(0,1);columnNorms=zeros(0,1);
else
    R=X*Z-Z*diag(returnedScaled);
    globalResidual=norm(R,'fro')/max(realmin, ...
        xnorm*norm(Z,'fro')+norm(Z*diag(returnedScaled),'fro'));
    Y=[Z,rightpairs_anti(Z)];singular=svd(Y);
    rankTolerance=max(size(Y))*eps*singular(1);
    basisRank=sum(singular>rankTolerance);
    rcondBasis=singular(end)/singular(1);
    columnNorms=sqrt(sum(abs(Z).^2,1)).';
end
pairingValid=all(pairResidual<=pairTolerance);
spectrumValid=pairingValid&&outputFinite&&~outputUnderflow;
eigenpairsValid=spectrumValid&&all(isfinite(colResidual))&& ...
    all(colResidual<=residualTolerance)&&isfinite(globalResidual)&& ...
    globalResidual<=residualTolerance;
fullBasis=basisRank==2*n;
rankH=NaN;if mod(basisRank,2)==0,rankH=basisRank/2;end
status='ok';message='Right eigenpairs checked; a numerically full quaternionic basis was observed.';
if ~eigenpairsValid
    status='undetermined';message='Conjugate pairing, representability or eigenpair residual checks failed.';
elseif ~fullBasis
    status='rankDeficient';
    message=['Right eigenpairs have small residuals, but V has no numerically full basis. ', ...
        'The input may be defective or too ill-conditioned; this is not a proof of defectiveness.'];
end
info=struct('schema','quaternion.right-eigenpairs.1','status',status, ...
    'isValid',eigenpairsValid,'message',message,'spectrumValid',spectrumValid, ...
    'eigenpairsValid',eigenpairsValid,'fullBasisObserved',fullBasis, ...
    'quaternionOrder',n,'algebraicCount',n,'matrixScale',scale, ...
    'scaledFullSpectrum',mu,'scaledRepresentatives',lambdaScaled,'returnedScaledRepresentatives',returnedScaled, ...
    'pairedIndices',pairs,'pairingResidual',pairResidual, ...
    'pairingTolerance',pairTolerance,'pairingValid',pairingValid, ...
    'selectedComplexColumns',selected,'selectedColumnsReflected',reflected, ...
    'clusterIndex',clusterIndex,'clusterTolerance',clusterTolerance, ...
    'selectionNovelty',novelty,'selectionTolerance',selectionTolerance, ...
    'columnResiduals',colResidual,'eigenResidual',globalResidual, ...
    'residualTolerance',residualTolerance,'columnNorms',columnNorms, ...
    'complexBasisRank',basisRank,'quaternionBasisRank',rankH, ...
    'basisSingularValues',singular,'rankTolerance',rankTolerance, ...
    'basisReciprocalCondition',rcondBasis,'outputFinite',outputFinite,'outputUnderflow',outputUnderflow, ...
    'spectrumType','right','normalization','unit Euclidean norm of each quaternion column', ...
    'assurance','floating-point; residual/rank checks are not interval proofs');
end

function y=rightpairs_anti(x)
% Anti-linear symplectic partner: eigenvalue mu -> conjugate(mu).
n=size(x,1)/2;y=[-conj(x(n+1:end,:));conj(x(1:n,:))];
end
