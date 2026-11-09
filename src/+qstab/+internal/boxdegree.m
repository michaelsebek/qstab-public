function [radius,witness] = boxdegree(M,variable)
%BOXDEGREE Exact componentwise regularity event (in exact arithmetic).
idx=M.grade+1;if strcmp(variable,'zinv'),idx=1;end
q=M.C(:,idx);w=M.weights(:,idx);radius=Inf;witness=[];
if any(w==0 & q~=0),return;end
ratios=zeros(4,1);mask=w>0;ratios(mask)=abs(q(mask))./w(mask);radius=max(ratios);
if ~isfinite(radius),return;end
C=M.C;C(:,idx)=0;D=C-M.C;
witness=struct('kind','degreeLoss','valid',true,'variable',variable,'grade',M.grade, ...
    'leadingIndexInInput',idx,'coefficientComponents',C,'perturbationComponents',D, ...
    'polynomial',struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:)), ...
    'weightedNorm',radius,'weights',M.weights,'evaluationResidual',NaN, ...
    'boundaryComponents',[],'validated',false);
end
