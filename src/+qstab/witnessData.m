function [C,Info] = witnessData(data,R,varargin)
%WITNESSDATA Reconstruct and independently check a returned family member.
%   [C,Info] = qstab.witnessData(data,R,...)
%   C holds the reconstructed 4-by-(grade+1) component array, or [] when
%   unavailable. Inspect Info.isAdmissible and Info.checkStatus before use.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,W] = qstab.witnessData(ORIGINAL_DATA,R) retains shared parameters/weights and
% checks coefficient membership before requesting a nominal spectral check.
% No optimization or robust solver is run. The original R is never modified.
% A degree-loss witness is NOT labelled an unstable finite root.
% Input/test variables and declared grade are explicit. In zinv, the nominal
% spectral check uses forward z, whereas inputRepresentatives are reciprocals
% in the original delay variable, with zero mapped to infinity (no clipping).
% A near-boundary spectrum may remain undetermined. Numerical membership is
% checked to tolerance, not claimed as a directed-rounding certificate.
%
% Options: Uncertainty='' (default from R), Tolerance=1e-10.
%
%   Examples (ascending coefficients):
%   F = {[-1 1],[-2 1]}; [~,R] = qstab.robuststable(F);
%   [C,W] = qstab.witnessData(F,R);
%   disp(C); disp(W.checkStatus);
%
%   Interpretation:
%   Reconstruction uses the original family and its shared parameters. A degree-loss witness is not a finite unstable zero. Membership and spectral checks are numerical, not outward-rounded proofs.
%
%   Reference: docs/reference/witnessData.html.
%   See also qstab.about.
try
if qstab.internal.isFamilyDefinition(data)
    if ~isstruct(R)||~isfield(R,'variable')||~isequal(data.variable,R.variable)
        error('qstab:family:ReportConflict','The saved result and the family must use the same variable.');
    end
    if isfield(R,'degree') && ~isequal(R.degree,data.grade)
        error('qstab:family:ReportConflict','The saved result and the family must use the same grade.');
    end
end
Info = qstab.internal.witnessDataReport(data,R,varargin{:});
C = Info.components;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
