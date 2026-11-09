function [A,info] = cubefamily(F)
%CUBEFAMILY qB02 eight-parameter Kharitonov family, exact before rounding.
L=F.intervalLower;U=F.intervalUpper;n=size(L,2)-1;
B=zeros(size(L));G=zeros(4,n+1,8);
parts=struct('evenMinus',zeros(size(L)),'evenPlus',zeros(size(L)), ...
    'oddMinus',zeros(size(L)),'oddPlus',zeros(size(L)));
for k=0:n
    parity=mod(k,2);sgn=(-1)^floor(k/2);
    if sgn>0,lo=L(:,k+1);hi=U(:,k+1);else,lo=U(:,k+1);hi=L(:,k+1);end
    B(:,k+1)=lo;
    for nu=1:4,G(nu,k+1,nu+4*parity)=hi(nu)-lo(nu);end
    if parity==0
        parts.evenMinus(:,k+1)=lo;parts.evenPlus(:,k+1)=hi;
    else
        parts.oddMinus(:,k+1)=lo;parts.oddPlus(:,k+1)=hi;
    end
end
A=F;A.kind='affine';A.base=B;A.generators=G;A.bounds=repmat([0,1],8,1);
A.vertices=[];A.intervalLower=[];A.intervalUpper=[];
A.parameterIndex=[(1:4).',zeros(4,1);(1:4).',ones(4,1)];
A.parameterNames={'e_real';'e_i';'e_j';'e_k';'o_real';'o_i';'o_j';'o_k'};
active=[];
for j=1:8,if any(any(G(:,:,j)~=0)),active(end+1)=j;end;end %#ok<AGROW>
leadNearest=min(max(zeros(4,1),L(:,end)),U(:,end));
info=struct('source','qB02 Kharitonov cube theorem','variable','s', ...
    'degree',n,'parameterCount',8,'activeParameters',active, ...
    'activeParameterCount',numel(active),'originalParameterCount',size(F.bounds,1), ...
    'parameterComponents',[1:4,1:4],'parameterParity',[zeros(1,4),ones(1,4)], ...
    'parameterNames',{A.parameterNames},'parts',parts, ...
    'sameFixedFrequencyValueSets',true,'subsetOfOriginalInterval',true, ...
    'degreeInvariant',any(leadNearest~=0),'leadingNormLowerBound',norm(leadNearest), ...
    'vertexOnlyTest',false,'representable',all(isfinite(B(:)))&&all(isfinite(G(:))), ...
    'validated',false);
end
