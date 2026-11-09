function R = leadbound(F,variable,tol,maxIterations,maxTime)
%LEADBOUND Supporting-hyperplane bound for the FULL leading-coefficient set.
% The separation proof uses an explicit support minimum, NOT convergence of
% the minimum-norm iteration. A positive bound therefore remains sufficient
% even when the Frank-Wolfe iterate is not an optimum. A near-zero residual
% is never promoted to an exact degree-loss witness.
started=tic;
idx=F.degree+1; testVariable=variable;
if strcmp(variable,'zinv'), idx=1; testVariable='z'; end
R=struct('status','undetermined','isSeparated',false,'lowerBound',0, ...
    'upperBound',Inf,'normalizedLowerBound',0,'scale',1,'direction',[], ...
    'iterations',0,'witness',[],'message','','coefficientIndex',idx-1, ...
    'variable',variable,'testVariable',testVariable,'validated',false, ...
    'assurance','floatingPointSupportingHyperplane','elapsedTime',0);
if strcmp(F.kind,'interval')
    L=F.intervalLower(:,idx); U=F.intervalUpper(:,idx);
    q=min(max(zeros(4,1),L),U); d=norm(q);
    if all(q==0)
        C=.5*F.intervalLower+.5*F.intervalUpper; C(:,idx)=q;
        theta=zeros(size(F.bounds,1),1);
        for j=1:numel(theta)
            theta(j)=C(F.parameterIndex(j,1),F.parameterIndex(j,2)+1);
        end
        R.status='zeroWitness'; R.witness=localWitness(C,[],theta,variable);
        R.message='The full Cartesian leading box contains zero.';
    elseif isfinite(d) && d>0
        R.status='separated';R.isSeparated=true;R.lowerBound=d;R.upperBound=d;
        R.scale=max(abs([L;U]));R.normalizedLowerBound=d/R.scale;
        R.direction=q/d;R.message='Exact Cartesian-box distance formula (floating-point evaluation).';
    else
        R.message='Cartesian-box distance is outside the representable range.';
    end
    R.elapsedTime=toc(started);return
end
if strcmp(F.kind,'polytope')
    A=reshape(F.vertices(:,idx,:),4,[]); m=size(A,2);
    weights=ones(m,1)/m; theta=[];
    if all(all(A==repmat(A(:,1),1,m)))
        q=A(:,1); d=norm(q);
        if all(q==0)
            weights=zeros(m,1);weights(1)=1;
            R.status='zeroWitness';
            R.witness=localWitness(F.vertices(:,:,1),weights,[],variable);
            R.message='The common leading coefficient is zero.';
        elseif isfinite(d) && d>0
            R.status='separated';R.isSeparated=true;R.lowerBound=d;R.upperBound=d;
            R.scale=max(abs(q));R.normalizedLowerBound=d/R.scale;
            R.direction=q/d;R.message='Common nonzero leading quaternion.';
        else
            R.message='Common leading norm is outside the representable range.';
        end
        R.elapsedTime=toc(started);return
    end
    scale=max(abs(A(:)));
    if ~(scale>0) || ~isfinite(scale)
        R.message='Leading coefficient scaling is unavailable.';return
    end
    As=A/scale;
else
    weights=[]; B=F.bounds;
    theta=.5*B(:,1)+.5*B(:,2);
    radius=.5*B(:,2)-.5*B(:,1);
    G=reshape(F.generators(:,idx,:),4,[]);
    c=F.base(:,idx)+G*theta;
    H=bsxfun(@times,G,radius.');
    if any(~isfinite(c)) || any(~isfinite(H(:)))
        R.message='Leading affine enclosure overflowed.';return
    end
    scale=max([abs(c);sum(abs(H),2)]);
    if scale==0
        C=localAffine(F,theta);
        if all(C(:,idx)==0)
            R.status='zeroWitness';R.witness=localWitness(C,[],theta,variable);
            R.message='The affine family has identically zero leading coefficient.';
        else
            R.message='Inconsistent near-zero affine leading computation.';
        end
        R.elapsedTime=toc(started);return
    end
    if ~isfinite(scale)
        R.message='Leading affine scaling is unavailable.';return
    end
    Cs=c/scale; Hs=H/scale; eta=zeros(numel(theta),1);
end
R.scale=scale;
for step=0:maxIterations
    if strcmp(F.kind,'polytope')
        q=As*weights; C=localConvex(F.vertices,weights);
    else
        q=Cs+Hs*eta; theta=.5*B(:,1)+.5*B(:,2)+radius.*eta;
        theta=min(max(theta,B(:,1)),B(:,2));
        C=localAffine(F,theta);
    end
    R.iterations=step;
    % Require exact zero of the reconstructed INPUT coefficient, not merely
    % a small normalized residual or a projected approximate feasible point.
    if all(C(:,idx)==0)
        R.status='zeroWitness';R.witness=localWitness(C,weights,theta,variable);
        R.message='An admissible convex/affine combination has a zero leading quaternion.';
        break
    end
    nq=norm(q);R.upperBound=nq*scale;
    if ~(nq>0) || ~isfinite(nq)
        R.message='No resolved direction for leading-coefficient separation.';break
    end
    u=q/nq;
    if strcmp(F.kind,'polytope')
        supports=u.'*As;[support,j]=min(supports); v=As(:,j);
        magnitude=max(sqrt(sum(As.^2,1)));
    else
        supports=Hs.'*u; etaVertex=-sign(supports);
        support=u.'*Cs-sum(abs(supports));v=Cs+Hs*etaVertex;
        magnitude=norm(Cs)+sum(sqrt(sum(Hs.^2,1)));
    end
    pad=(tol+128*eps)*max(1,magnitude);
    lb=support-pad;
    if isfinite(lb) && lb>0 && isfinite(lb*scale) && lb*scale>0
        R.status='separated';R.isSeparated=true;R.lowerBound=lb*scale;
        R.normalizedLowerBound=lb;R.direction=u;
        R.message='A supporting hyperplane strictly separates the full leading set from zero.';
        break
    end
    if step>=maxIterations
        R.message='Leading-coefficient iteration budget reached without a separator.';break
    end
    if toc(started)>=maxTime
        R.message='Time budget reached during leading-coefficient separation.';break
    end
    d=v-q;dd=d.'*d;
    if ~(dd>0) || ~isfinite(dd)
        R.message='Leading-coefficient iteration stagnated without a strict separator.';break
    end
    gamma=min(1,max(0,-(q.'*d)/dd));
    if ~(gamma>0)
        R.message='Leading-coefficient iteration has no resolved descent step.';break
    end
    if strcmp(F.kind,'polytope')
        weights=(1-gamma)*weights;weights(j)=weights(j)+gamma;
        weights=weights/sum(weights);
    else
        eta=(1-gamma)*eta+gamma*etaVertex;
        eta=min(1,max(-1,eta));
    end
end
R.elapsedTime=toc(started);
end

function C=localAffine(F,t)
C=F.base;
for j=1:numel(t),C=C+t(j)*F.generators(:,:,j);end
end
function C=localConvex(V,w)
C=zeros(size(V,1),size(V,2));
for j=1:numel(w),C=C+w(j)*V(:,:,j);end
end
function W=localWitness(C,w,t,variable)
Ct=C;tv=variable;
if strcmp(variable,'zinv'),Ct=C(:,end:-1:1);tv='z';end
W=struct('kind','zeroLeadingCoefficient','components',Ct, ...
    'originalComponents',C,'originalWeights',w,'parameter',t, ...
    'coefficientVariable',tv,'originalCoefficientVariable',variable);
end
