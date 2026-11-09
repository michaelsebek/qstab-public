function [membership,Info] = ballset(P,weights,z,r,varargin)
%BALLSET Pointwise spectral-set membership for a weighted coefficient ball.
%   [MEMBERSHIP,Info] = qstab.ballset(P,weights,z,r,...)
%   MEMBERSHIP is pointwise: 1 inside, 0 outside, NaN unresolved.
%   It is not a whole-family stability decision.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,S] = qstab.ballset(P,GAMMA,Z,R,'Exponent',ETA) classifies the requested complex
% points only. Inside means some admissible member has this algebraic class;
% outside means the pointwise distance exceeds R with a numerical guard.
% Values in the tolerance band remain unresolved. No grid certifies the
% continuum between its points or robust stability of the whole family.
% Z refers to the INPUT polynomial variable, including a delay variable.
% Options are passed to qstab.ballvalue. Witness generation is optional.
%
%   Examples (ascending coefficients):
%   [membership,S] = qstab.ballset([2 1],[1 0],[-2 -1 0],0.5);
%   disp(membership);
%
%   Interpretation:
%   Only the supplied query points are classified. A sampled contour or grid never certifies the whole stability region.
%
%   Reference: docs/reference/ballset.html.
%   See also qstab.about.
try
Info = qstab.internal.ballsetReport(P,weights,z,r,varargin{:});
membership = NaN(size(Info.inside));
membership(Info.inside) = 1; membership(Info.outside) = 0;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
