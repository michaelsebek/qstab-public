function [points,Info] = valueset(P, lambda, varargin)
%VALUESET Evaluate scalar quaternion norm-polynomial value sets.
%   [POINTS,Info] = qstab.valueset(P, lambda, ...)
%   POINTS is the complex array of norm values (samples by arguments).
%   Sampling alone does not prove zero exclusion.
%   Info is the full diagnostic structure returned by earlier releases.
%   [~,S] = qstab.valueset(P,LAMBDA) evaluates N_P at the finite complex points
%   LAMBDA. P is a scalar polynomial in the ascending coefficient format
%   accepted by qstab.polycomponents. If P = p0+p1*i+p2*j+p3*k, where the
%   p_j are REAL coefficient polynomials in a central variable, then
%
%       N_P(lambda) = p0(lambda)^2 + p1(lambda)^2 ... + p3(lambda)^2.
%
%   These are ordinary complex squares, NOT squared absolute values.
%   In particular, this is not abs(P(lambda))^2 or a quaternion point
%   evaluation at a chosen imaginary direction.
%
%   [~,S] = qstab.valueset({P1,...,Pm},LAMBDA,'Weights',W) first forms the
%   coefficientwise convex combinations sum_j W(k,j)*Pj. A row of the real
%   matrix W is one parameter sample. The SAME row is used for every
%   coefficient and quaternion component. The norm map is quadratic, so
%   these values are not convex combinations of the vertex norm values.
%   Without W, only the m vertices are evaluated (W = eye(m)).
%   A cell of scalar coefficients WITHOUT the Weights option is treated as
%   one polynomial, consistently with qstab.polycomponents. For a family
%   of constant polynomials, specify Weights explicitly (possibly []).
%
%   Output S.points(k,j) corresponds to parameter row k and LAMBDA(j).
%   S.parameters contains the weights actually used. Weights within a
%   floating-point tolerance of the simplex are clipped/normalized onto
%   it; invalid weights raise an error. S.inputWeights records the input.
%   S.kind is 'singleton' for a fixed polynomial and 'samples' for a
%   nontrivial vertex family. S.isCertificate is always false: a finite
%   cloud of values does not prove zero exclusion over a continuum.
%
%   Examples:
%       [~,S] = qstab.valueset([1 1],1i); % P(s)=1+s, N_P(1i)=2i
%       t = linspace(0,1,21).';
%       [~,S] = qstab.valueset({[1 2 1],[2 3 1]},1i*[0 .5 1], ...
%                         'Weights',[1-t t]);
%       qstab.plotValueSet(S);
%
%   This function applies to SCALAR polynomials. It does not use an
%   entrywise norm as a substitute for a polynomial-matrix determinant.
%   See also qstab.polycomponents, qstab.plotValueSet.
%
%   Examples (ascending coefficients):
%   [points,S] = qstab.valueset({[1 1],[2 1]},1i,'Weights',[1 0;.5 .5;0 1]);
%   disp(points);
%
%   Interpretation:
%   Weights combine all original polynomial coefficients together before evaluating the quadratic real norm map. These are attainable samples, not an outer enclosure or a zero-exclusion proof.
%
%   Reference: docs/reference/valueset.html.
%   See also qstab.about.
try
if qstab.internal.isFamilyDefinition(P)
    [P,kind]=qstab.internal.familyDefinition(P);
    if ~strcmp(kind,'polytope'),error('qstab:family:ValuesetKind','valueset accepts a polynomial or a polytope; for other families use plotFamilyValues with Plot=false.');end
    if isempty(varargin),varargin={'Weights',[]};end
end
Info = qstab.internal.valuesetReport(P, lambda, varargin{:});
points = Info.points;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
