function [window,nominalAdmissible,status]=parameterRay(B,u)
%PARAMETERRAY Closed, geometrically admissible parameter interval on t=alpha*u.
% This does not detect any stability contact. B=[] means the entire H.
if isempty(B),window=[-Inf Inf];nominalAdmissible=true;status='unrestricted';return;end
u=double(u(:));B=double(B);
nominalAdmissible=all(B(:,1)<=0 & B(:,2)>=0);
if any((u==0)&(B(:,1)>0|B(:,2)<0))
    window=[];status='empty';return
end
active=u~=0;rows=B(active,:);den=u(active);R=bsxfun(@rdivide,rows,den);
if any(~isfinite(R(:)))||any(rows(:)~=0&R(:)==0)
    window=[NaN NaN];status='notRepresentable';return
end
lo=max(min(R,[],2));hi=min(max(R,[],2));
if lo>hi,window=[];status='empty';else,window=[lo hi];status='closed';end
end
