function [regions,Info] = certificateEnclosure(data,R,varargin)
%CERTIFICATEENCLOSURE Reconstruct ONE retained value-space certificate cell.
%   [REGIONS,Info] = qstab.certificateEnclosure(data,R,...)
%   REGIONS is the cell array of local outer regions, with polygon and controls.
%   Use their union at this argument; Info records the scope and coverage.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,E] = qstab.certificateEnclosure(F,R,'CellIndex',j) reads R.certificates{j}.
% For a polytope, F must be the ORIGINAL family: original vertex weights are
% reused, including four-generator reductions. Delay coefficients are reversed
% at the recorded grade before evaluating the TEST variable z. This operation
% is explicit in E.originalVariable / E.variable.
% Affine and even/odd interval leaves already store rectangles; F may then be
% [] because their geometry is read from the report. Their source data are not
% revalidated. Lyapunov certificates have no norm-value image and are rejected.
%
% Bernstein controls or affine centre rectangles are inflated by their saved
% whole-frequency variation radius. A CIRCUMSCRIBED disk polygon is used for
% the display, so it is not a falsely shrunken inscribed approximation.
% This reads/reconstructs geometry; it never reruns the stability solver and
% never elevates a floating-point report to a validated interval proof.
%
%   Examples (ascending coefficients):
%   V = {[1 1],[2 1]};
%   [~,R] = qstab.zeroexclude(V,'StoreCertificates',true);
%   [regions,E] = qstab.certificateEnclosure(V,R,'CellIndex',1);
%
%   Interpretation:
%   Reconstructs an outer set for an entire retained parameter/frequency cell. It checks compatibility of the original family. Plotting a cell does not add global coverage.
%
%   Reference: docs/reference/certificateEnclosure.html.
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
Info = qstab.internal.certificateEnclosureReport(data,R,varargin{:});
regions = Info.regions;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
