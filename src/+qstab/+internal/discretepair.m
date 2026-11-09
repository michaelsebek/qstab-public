function M = discretepair(P0,P1,variable,grade)
%DISCRETEPAIR Preserve one quaternion parameter and declared temporal grade.
% No memberwise monic normalization; C0/C1 are the original physical data.
variable=lower(char(variable));
if ~any(strcmp(variable,{'z','zinv'}))
    error('qstab:discretepair:Variable','The discrete variable must be z or zinv.');
end
C0=qstab.polycomponents(P0);C1=qstab.polycomponents(P1);
g=max(size(C0,2),size(C1,2))-1;
if ~isempty(grade)
    if ~(isnumeric(grade)&&isreal(grade)&&isscalar(grade)&&isfinite(grade)&&grade==fix(grade)&&grade>=g)
        error('qstab:discretepair:Grade','Grade must be an integer at least as large as both supplied layouts.');
    end
    g=double(grade);
end
C0(:,end+1:g+1)=0;C1(:,end+1:g+1)=0;
M=struct('variable',variable,'testVariable','z','grade',g, ...
    'inputC0',C0,'inputC1',C1,'C0',C0,'C1',C1, ...
    'scale0',max(abs(C0(:))),'scale1',max(abs(C1(:))), ...
    'logScaleRatio',NaN,'S0',C0,'S1',C1,'hasPerturbation',any(C1(:)~=0));
if strcmp(variable,'zinv'),M.C0=fliplr(C0);M.C1=fliplr(C1);end
if M.scale0>0,M.S0=M.C0/M.scale0;else,M.S0=M.C0;end
if M.scale1>0,M.S1=M.C1/M.scale1;else,M.S1=M.C1;end
if any(M.C0(:)~=0 & M.S0(:)==0)||any(M.C1(:)~=0 & M.S1(:)==0)
    error('qstab:discretepair:DynamicRange','Scaling would erase a nonzero coefficient; rescale the supplied data.');
end
if M.scale0>0&&M.scale1>0,M.logScaleRatio=log(M.scale0)-log(M.scale1);end
end
