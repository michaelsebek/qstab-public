function [W,D] = boxcandidate(M,C,r,variable,tol)
%BOXCANDIDATE Negative-only check of one member of the ORIGINAL component box.
% A computed or padded endpoint can be just outside the requested box. Move
% adjustable components inward, then check BOTH abs(delta)<=r*weights and
% max(abs(delta)./weights)<=r, with NO membership-tolerance enlargement.
% Fixed components must already agree exactly; they are never silently fixed.
% An inward move does not preserve instability: recheck the resulting member.
% Stable sampled members are NOT certificates for the complete family.
W=[];
D=struct('schema','qstab.box-candidate.1','status','undetermined', ...
    'insideRequestedBox',false,'fixedComponentsPreserved',false, ...
    'retracted',false,'inwardSteps',0,'fallbackToCenter',false, ...
    'originalComponents',C,'coefficientComponents',[], ...
    'weightedNorm',Inf,'requestedRadius',r,'nominal',[], ...
    'spectralCheck',[],'negativeScreeningOnly',true,'message','');
if ~isnumeric(C)||~isreal(C)||~isequal(size(C),size(M.C))||any(~isfinite(C(:)))
    D.message='Candidate coefficients have an invalid shape or nonfinite value.';return
end
if ~(isnumeric(r)&&isreal(r)&&isscalar(r)&&isfinite(r)&&r>=0)
    D.message='Candidate radius must be finite and nonnegative.';return
end
active=M.weights>0;
if any(C(~active)~=M.C(~active))
    D.message='A candidate changed a fixed component; it is rejected, not repaired.';return
end
D.fixedComponentsPreserved=true;
span=r*M.weights;
loEndpoint=M.C-span;hiEndpoint=M.C+span;
if any(~isfinite([span(:);loEndpoint(:);hiEndpoint(:)])) || ...
        (r>0 && any(active(:) & span(:)==0))
    D.message='Original box endpoints or positive widths are not representable.';return
end
original=C;
C(active)=min(max(C(active),loEndpoint(active)),hiEndpoint(active));
% Endpoint arithmetic can round outward. Walk only offending coordinates
% toward their exact stored centre. eps(abs(x)) is a conservative step at
% a binade transition; skipping one representable neighbour is harmless.
for pass=1:8
    [inside,weighted,bad]=localMembership(M,C,span,r,active);
    if inside,break;end
    indices=find(bad);
    for k=1:numel(indices)
        j=indices(k);step=eps(abs(C(j)));
        if C(j)>M.C(j)
            C(j)=max(M.C(j),C(j)-step);
        elseif C(j)<M.C(j)
            C(j)=min(M.C(j),C(j)+step);
        else
            C(j)=M.C(j);
        end
    end
    D.inwardSteps=D.inwardSteps+numel(indices);
end
[inside,weighted,bad]=localMembership(M,C,span,r,active);
if ~inside
    % No tolerance relaxation: the stored centre is always admissible. Any
    % such extra movement must still survive the independent nominal check.
    C(bad)=M.C(bad);D.fallbackToCenter=true;
    [inside,weighted]=localMembership(M,C,span,r,active);
end
D.retracted=any(C(:)~=original(:));D.coefficientComponents=C;
D.insideRequestedBox=inside;D.weightedNorm=weighted;
if ~inside
    D.message='No finite candidate passed the original component bounds.';return
end
lead=M.grade+1;if strcmp(variable,'zinv'),lead=1;end
isDegreeLoss=all(C(:,lead)==0);
if isDegreeLoss
    proposed='assumptionFailed';
else
    D.nominal=qstab.internal.isstableReport(localParts(C),'Variable',variable, ...
        'Grade',M.grade,'Tolerance',tol);
    if strcmp(D.nominal.status,'stable')
        D.status='stable';
        D.message='This one candidate is stable; no complete-box conclusion follows.';return
    end
    proposed='notStable';
end
Q=struct('status',proposed,'variable',variable, ...
    'method','originalComponentBoxCandidate', ...
    'witness',struct('originalComponents',C));
% witnessData performs an INDEPENDENT nominal spectral test. Its generic
% family-membership tolerance is not used to enlarge this original box:
% the exact stored components already passed the strict checks above.
S=qstab.internal.witnessDataReport({localParts(loEndpoint),localParts(hiEndpoint)},Q, ...
    'Uncertainty','interval','Tolerance',tol);
D.spectralCheck=S;
if strcmp(proposed,'notStable') && S.isConfirmedCounterexample
    kind='unstableMember';D.status='notStable';
elseif strcmp(proposed,'assumptionFailed') && S.isAdmissible && S.isDegreeLoss
    kind='degreeLoss';D.status='assumptionFailed';
else
    D.message='An admissible candidate has no independently confirmed negative decision.';return
end
delta=C-M.C;
W=struct('kind',kind,'valid',true,'variable',variable,'grade',M.grade, ...
    'coefficientComponents',C,'perturbationComponents',delta, ...
    'polynomial',localParts(C),'weightedNorm',weighted,'weights',M.weights, ...
    'spectralCheck',S,'requestedRadius',r,'insideRequestedBox',true, ...
    'fixedComponentsPreserved',true,'boundaryComponents',[], ...
    'retractedToOriginalBox',D.retracted,'validated',false);
D.message='A member inside the original requested box has an independently confirmed negative decision.';
end

function [inside,weighted,bad]=localMembership(M,C,span,r,active)
delta=C-M.C;ratios=zeros(size(C));
ratios(active)=abs(delta(active))./M.weights(active);
weighted=max(ratios(:));
bad=active & (~isfinite(delta) | ~isfinite(ratios) | abs(delta)>span | ratios>r);
inside=~any(bad(:)) && all(delta(~active)==0) && isfinite(weighted);
end

function P=localParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
