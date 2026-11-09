function [rows,Info] = certificateData(R)
%CERTIFICATEDATA Read retained robust-test cells, without solving again.
%   [ROWS,Info] = qstab.certificateData(R)
%   ROWS contains retained cell descriptors. Info contains counts and coverage.
%   Info is the full diagnostic structure returned by earlier releases.
% [~,D] = qstab.certificateData(R) accepts qstab.robuststable / zeroexclude reports.
% All accepted and unresolved records are kept, including compact unvisited
% subfamily descriptors. The reported solver verdict is NEVER recomputed.
% reportedComplete is the solver's claim; storageComplete concerns retained
% accepted leaves ONLY, not an independent proof of parameter coverage.
% An empty certificates list does not undo a solver's stable verdict, nor
% does a plotted frequency projection prove joint-parameter coverage.
% Metrics are labelled by their units; Lyapunov slack is NOT a spectral or
% coefficient-space robustness margin. No eig, roots or solvers are called.
%
%   Examples (ascending coefficients):
%   [~,R] = qstab.zeroexclude({[1 1],[2 1]},'StoreCertificates',true);
%   [rows,D] = qstab.certificateData(R);
%   disp(rows);
%
%   Interpretation:
%   This reads retained evidence and never reruns the robust solver. Absent stored leaves are not fabricated. Lyapunov slack, norm-value separation and parameter radius have different units and meanings.
%
%   Reference: docs/reference/certificateData.html.
%   See also qstab.about.
try
Info = qstab.internal.certificateDataReport(R);
rows = Info.rows;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
