function [ax,D] = plotFamilyValues(F,omega,varargin)
%PLOTFAMILYVALUES Plot values of sampled family members at one frequency.
%   qstab.plotFamilyValues(F,omega)                100 sampled members
%   qstab.plotFamilyValues(F,omega,N)              N sampled members
%   [ax,D] = qstab.plotFamilyValues(F,omega,N,...) plot and keep the data
%   [~,D] = qstab.plotFamilyValues(F,omega,N,'Plot',false)  no figure
%
%   For every sampled P(v)=sum a_k v^k (coefficients on the left), compute
%   its RIGHT point evaluation w=sum a_k*xi^k. Choose a pure unit quaternion
%   u with 'AxisUnit' (default i). For variable s use xi=u*omega, where
%   omega is in radians per unit time. For z use cos(omega)+u*sin(omega),
%   for delay use cos(omega)-u*sin(omega); omega is in radians PER SAMPLE
%   in both discrete cases. Signed real omega is accepted. To query physical
%   frequency Omega in rad/s at sample period h, pass omega=Omega*h.
%
%   'Projection'   'class' (default): (Re w,norm(Im w)). This preserves
%                  distance to the origin and maps only zero to the origin;
%                  directions in the imaginary 3-space are not preserved.
%                  'slice': (Re w,Im(w) dot u), colored by the length of
%                  its discarded imaginary component perpendicular to u.
%                  A point at the plotted origin may have w ~= 0.
%                  'normpoly': planar complex N_P(zeta)=sum_l h_l(zeta)^2
%                  for the FOUR real component polynomials h_l, at central
%                  complex zeta. This tests algebraic zero classes; it is
%                  not the quaternion value w or its squared norm.
%   'AxisUnit'     'i','j','k', a real unit 3-vector, or a scalar pure unit
%                  MATLAB quaternion or a,b,c,d component struct. Fixes the evaluation direction;
%                  the norm-polynomial plane is independent of this choice.
%   'AxisUnit','sphere' with 'Projection','class' samples K input directions
%                  per polynomial: +i,+j,+k,-i,-j,-k first (as many as fit),
%                  then a Fibonacci sphere. 'DirectionCount',K defaults to 32.
%                  N counts members; N*K <= 1e6. Values at real arguments
%                  repeat across all directions. No sphere search is needed
%                  in normpoly mode. A signed slice requires one fixed axis.
%   F may be a struct from qstab.family or qstab.randfamily, or any raw
%   input accepted by qstab.familycomponents.
%   By default, a cell {P1,...,Pm} denotes a polytope. To sample independent
%   component bounds {L,U}, pass 'Uncertainty','interval'. Affine and
%   oneparameter definitions are structs (see qstab.familycomponents).
%   A canonical qstab.familycomponents result is accepted as well.
%   'Variable'      's' (default), 'z' or 'zinv'; matching family metadata
%                   are inferred for tagged family structs.
%   'Grade'         For raw inputs, pad to this grade (never discard coefficients).
%   'Uncertainty'   'polytope', 'interval', 'affine' or 'oneparameter'.
%
%   'Sampling'     'grid' (default) or 'random'; N counts family members.
%   'Seed'         Integer 0..2^32-1 for repeatable random sampling with
%                  no change to the caller's RNG. Without Seed, random
%                  sampling advances the caller's RNG.
%   N is at most 1e6, with at most 1e7 coefficient/parameter samples.
%   'ParameterBox' Finite 4-by-2 component limits for a oneparameter model;
%                  required for an unrestricted t; a subbox for bounded t.
%   'Plot'         true (default) or false; false returns ax=[].
%   'Parent'       Existing axes; 'Visible','off' hides a new figure.
%
%   D.valueComponents has N*K rows (K=1 by default), with columns
%   [real i j k] of w. D.memberIndex and D.directionIndex identify each row.
%   D.valueNorm is the Euclidean norm of each w. D.projected has two
%   columns. D.offSlice records the
%   perpendicular imaginary length. D.normValues is the complex N-by-1
%   array in 'normpoly' mode, empty otherwise. D.complexArgument is zeta,
%   D.argumentComponents lists the K query quaternions xi; D.axisUnit lists
%   their K pure unit directions as real 3-vectors.
%   D.sample gives weights/parameters; D.coefficients(:,k+1,n) contains the
%   four real parts of member n at power k. D.isCertificate=false: these
%   are samples, not a certified value enclosure or zero-exclusion test.
%
%   Example (ascending coefficients; no add-ons needed):
%     minusJ = struct('a',[0 1],'b',[0 0],'c',[-1 0],'d',[0 0]);
%     F = qstab.family('polytope',{[0 1],minusJ});
%     qstab.plotFamilyValues(F,1,41,'AxisUnit','j');          % actual P(j)
%     qstab.plotFamilyValues(F,1,41,'Projection','normpoly'); % central N_P(1i)
%   See also qstab.plotFamilyZeros, qstab.valueset, qstab.plotValueSet, qstab.valueenclosure.
try
    if nargin<2 || ~(isnumeric(omega)&&isreal(omega)&&isscalar(omega)&&isfinite(omega))
        error('qstab:sample:Frequency','Supply a scalar polynomial family and one finite real angular frequency.');
    end
    [n,o] = qstab.internal.familyPlotOptions(varargin,true);
    [C,S] = qstab.internal.sampleFamily(F,n,o);
    u=o.axis; omega=double(omega);
    if strcmp(S.variable,'s')
        z=complex(0,omega);
    else
        z=complex(cos(omega),sin(omega));
        % The canonical real boundary angles 0 and +/-pi specify exact
        % points, irrespective of the floating approximation to pi.
        if omega==0, z=1; elseif abs(omega)==pi, z=-1; end
        if strcmp(S.variable,'zinv'), z=conj(z); end
    end
    % Evaluate four real component polynomials in central complex arithmetic.
    % If h_l(z)=A_l+i_c*B_l, the actual quaternion evaluation is A+B*u,
    % not A+u*B. The ordering matters for coefficients outside the u slice.
    H=reshape(C(:,end,:),4,n);
    for power=S.grade:-1:1
        H=H*z+reshape(C(:,power,:),4,n);
    end
    H=H.'; A=real(H); B=imag(H);
    ndir=size(u,1);
    memberIndex=kron((1:n).',ones(ndir,1));
    directionIndex=repmat((1:ndir).',n,1);
    directions=u(directionIndex,:);
    A=A(memberIndex,:); B=B(memberIndex,:);
    V=A(:,2:4)+bsxfun(@times,B(:,1),directions)+cross(B(:,2:4),directions,2);
    W=[A(:,1)-sum(B(:,2:4).*directions,2),V];
    rad=hypot(hypot(V(:,1),V(:,2)),V(:,3));
    signed=sum(V.*directions,2); normW=hypot(W(:,1),rad);
    perpendicular=cross(V,directions,2);
    rem=hypot(hypot(perpendicular(:,1),perpendicular(:,2)),perpendicular(:,3));
    if any(~isfinite([W(:);normW;rem;signed]))
        error('qstab:sample:ArithmeticRange','A quaternion evaluation is not representable. Rescale the polynomial or frequency.');
    end
    nv=[];
    switch o.projection
        case 'class', xy=[W(:,1),rad];
        case 'slice', xy=[W(:,1),signed];
        case 'normpoly'
            % Scale only when squaring could overflow or vanish. Keeping
            % ordinary Horner values for typical magnitudes avoids adding
            % roundoff around exact cancellations such as N_(s-j)(1i)=0.
            sc=max(max(abs(real(H)),abs(imag(H))),[],2);
            scaled=sc>sqrt(realmax/16) | (sc~=0 & sc<sqrt(realmin));
            nv=complex(zeros(n,1)); nv(~scaled)=sum(H(~scaled,:).^2,2);
            if any(scaled)
                hs=bsxfun(@rdivide,H(scaled,:),sc(scaled));
                v=sum(hs.^2,2); r=(v.*sc(scaled)).*sc(scaled); nv(scaled)=r;
                if any(v~=0 & r==0)
                    error('qstab:sample:ArithmeticRange','A nonzero norm-polynomial value underflowed. Rescale the polynomial.');
                end
            end
            if any(~isfinite(nv))
                error('qstab:sample:ArithmeticRange','A norm-polynomial value overflowed. Rescale the polynomial.');
            end
            xy=[real(nv),imag(nv)];
    end
    D=struct('sample',S,'coefficients',C,'frequency',omega,'axisUnit',u, ...
        'complexArgument',z,'argumentComponents',[repmat(real(z),ndir,1),imag(z)*u], ...
        'memberIndex',memberIndex,'directionIndex',directionIndex, ...
        'projection',o.projection,'valueComponents',W,'valueNorm',normW, ...
        'offSlice',rem,'normValues',nv,'projected',xy,'isCertificate',false);
    ax=[];
    if o.plot, ax=qstab.internal.drawSampledFamily('values',D,o); end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
