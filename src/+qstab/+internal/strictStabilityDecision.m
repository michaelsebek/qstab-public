function tf = strictStabilityDecision(Info)
%STRICTSTABILITYDECISION Single public-layer mapping of strict stability status.
% Numerical stability tests and private object adapters must not round or
% otherwise reinterpret this decision at their interfaces.
switch Info.status
    case 'stable', tf=1;
    case 'notStable', tf=0;
    case {'undetermined','assumptionFailed'}, tf=NaN;
    otherwise
        error('qstab:DecisionStatus','Unexpected strict stability status: %s',Info.status);
end
end
