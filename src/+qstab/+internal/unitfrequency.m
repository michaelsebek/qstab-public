function [rho,info] = unitfrequency(P0,P1,theta,o)
%UNITFREQUENCY Vector wrapper of full-class discrete one-parameter quotients.
if any(theta(:)>pi),error('qstab:freqmargin:Angle','Discrete angles must lie in [0,pi].');end
M=qstab.internal.discretepair(P0,P1,lower(char(o.Variable)),o.Grade);
if M.scale0==0||all(M.C0(:,end)==0)
    error('qstab:freqmargin:Nominal','The nominal forward leading coefficient must be nonzero at the declared grade.');
end
if o.CheckNominal
    nominal=qstab.internal.isstableReport(parts(M.C0),'Variable','z','Grade',M.grade,'Tolerance',o.Tolerance);
    if ~strcmp(nominal.status,'stable'),error('qstab:freqmargin:Nominal','P0 must be numerically verified Schur.');end
else,nominal=struct('status','notChecked');end
rho=NaN(size(theta));nw=numel(theta);side=lower(char(o.Side));
info=struct('side',side,'variable',M.variable,'testVariable','z','grade',M.grade, ...
    'parameterMeaning','angle of forward z, radians','assurance','floating-point','validated',false, ...
    'valid',false(size(theta)),'parameterComponents',NaN(4,nw), ...
    'boundaryComponents',NaN(4,nw),'forwardBoundaryComponents',NaN(4,nw), ...
    'evaluationResidual',NaN(size(theta)),'objectiveRadius',NaN(size(theta)), ...
    'centerLowerEstimate',zeros(size(theta)),'minNominalModulus',NaN(size(theta)), ...
    'witnesses',{cell(size(theta))},'message',{cell(size(theta))},'nominal',nominal);
for k=1:nw
    E=qstab.internal.unitquotient(M,double(theta(k)),side,o.Tolerance);
    rho(k)=E.radius;info.valid(k)=E.valid;info.parameterComponents(:,k)=E.parameterComponents;
    info.boundaryComponents(:,k)=E.boundaryComponents;info.forwardBoundaryComponents(:,k)=E.forwardBoundaryComponents;
    info.evaluationResidual(k)=E.evaluationResidual;info.objectiveRadius(k)=E.radius;
    info.centerLowerEstimate(k)=E.lower;info.minNominalModulus(k)=E.scaledMinNominal*M.scale0;
    info.witnesses{k}=E.witness;info.message{k}=E.message;
end
end
function P=parts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
