function R = ballstableReport(P,weights,r,varargin)
%BALLSTABLEREPORT Internal implementation returning the full report.
if ~(isnumeric(r)&&isreal(r)&&isscalar(r)&&isfinite(r)&&r>=0)
    error('qstab:ballstable:Radius','Radius must be a finite nonnegative real scalar.');
end
M=qstab.internal.ballradiusReport(P,weights,varargin{:});
R=struct('schema','qstab.ballstable.1','status','undetermined','isStable',[], ...
    'radius',double(r),'margin',M,'witness',[],'variable',M.variable, ...
    'validated',false,'message','The requested radius is not separated from the computed global bracket.');
if isfield(M.nominal,'status')&&strcmp(M.nominal.status,'notStable')
    R.status='notStable';R.isStable=false;
    R.message='The nominal member itself is not stable and belongs to every ball.';
    C=qstab.polycomponents(P,'Grade',M.grade);
    Q=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
        'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
    R.witness=struct('kind','nominal','polynomial',Q,'weightedNorm',0,'valid',true);return
end
if strcmp(M.status,'assumptionFailed')
    if strcmp(M.nominal.status,'assumptionFailed'),R.status='assumptionFailed';end
    R.message=M.message;return
end
if r==0&&strcmp(M.nominal.status,'stable')
    R.status='stable';R.isStable=true;R.message='Radius zero is the verified nominal member.';return
end
if strcmp(M.status,'undetermined'),R.message=M.message;return;end
% An explicit degree-loss member uses the required forward leading
% coefficient for delay inputs. This check precedes finite-root claims.
if isfinite(M.degreeLossRadius)&&r>=M.degreeLossRadius
    R.status='assumptionFailed';R.witness=M.degreeLossWitness;
    R.message='The ball contains a member with a lost required leading coefficient.';return
end
% A closed ball includes its boundary. At numerical equality, accept
% only an independently identified exact boundary member with a provably
% admissible distance in the ORIGINAL coefficient coordinates.
if ~isempty(M.witness) && M.witness.valid && strcmp(M.witness.kind,'finiteZero') && ...
        M.witness.weightedNorm<=r && r<=M.upper
    C=qstab.polycomponents(P,'Grade',M.grade); Q=M.witness.coefficientComponents;
    if localAxisMember(C,Q,M.witness.weights,r)
        Ct=Q; domain='continuous';
        if ~strcmp(M.variable,'s'),domain='discrete';end
        if strcmp(M.variable,'zinv'),Ct=fliplr(Ct);end
        B=qstab.internal.exactBoundary('polynomial',Ct,domain);
        if ~isempty(B)
            R.status='notStable';R.isStable=false;R.witness=M.witness;
            R.witness.exactBoundary=B;
            R.message='The closed ball contains an exactly checked boundary member.';return
        end
    end
end
if r<M.lower
    R.status='stable';R.isStable=true;
    R.message='The requested radius is strictly below the numerical global lower bound.';return
end
if isfinite(M.upper)&&r>M.upper&&~isempty(M.witness)&&M.witness.valid&& ...
        strcmp(M.witness.kind,'finiteZero')&&M.witness.weightedNorm<=r
    R.status='notStable';R.isStable=false;R.witness=M.witness;
    R.message='An admissible numerically checked boundary-zero member lies inside the requested ball.';
end
end

function yes=localAxisMember(C,Q,weights,r)
% For any l_eta norm a perturbation of ONE real component has distance
% |q-c| / gamma_k. Use an error-free subtraction and a conservative
% product bound; an unproved equality must remain unresolved.
yes=false;
if ~isequal(size(C),size(Q)) || any(~isfinite(Q(:))),return;end
indices=find(C~=Q);if numel(indices)~=1,return;end
idx=indices(1); [~,j]=ind2sub(size(C),idx);w=weights(j);if w==0 || r==0,return;end
a=Q(idx);b=-C(idx);if abs(a)<abs(b),v=a;a=b;b=v;end
d=a+b;if ~isfinite(d) || b-(d-a)~=0,return;end
% Compare d with r*w using mantissas (their product cannot over/underflow).
[fd,ed]=log2(abs(d));[fr,er]=log2(r);[fw,ew]=log2(w);
m=fr*fw; % in [1/4,1)
% Multiplying a mantissa by 1/2 is exact; otherwise round DOWN by a
% full ulp, sufficient in round-to-nearest arithmetic.
if fr~=0.5 && fw~=0.5,m=m-eps(m);end
[fm,em]=log2(m);ep=er+ew+em;
yes=ed<ep || (ed==ep && fd<=fm);
end
