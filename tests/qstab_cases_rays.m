function cases = qstab_cases_rays()
%QSTAB_CASES_RAYS Analytic fixed-ray margins and numerical ambiguity tests.
%   Returns {description,@zeroArgumentTest}. These tests require only native
%   numeric/component data. Octave execution is auxiliary to native MATLAB.
    cases = { ...
        'Lyapunov operator has the analytic scalar and identity actions', @localOperatorIdentity; ...
        'qB01 cubic re-entry: both first contacts are selected', @localReentry; ...
        'qB01 noncommuting sides: finite first contacts survive infinite QZ modes', @localSides; ...
        'Zero and nilpotent directions give the correct infinite endpoints', @localInfinity; ...
        'A tiny nonzero direction has a distant finite first contact', @localTinyDirection; ...
        'Near-real complex contacts are not projected by the user tolerance', @localNearReal; ...
        'A defective touching contact cannot be skipped in favor of a later loss', @localTangency; ...
        'Nonmonic and extreme common coefficient scaling preserve the margins', @localScaling; ...
        'The direction is not normalized and parameter norm is reported separately', @localDirection; ...
        'Fixed-ray assumptions are checked before reporting a margin', @localAssumptions ...
    };
end

function localOperatorIdentity()
    scalar=localParts(-3,2,-5,7);
    [L,info]=qstab.lyapoperator(scalar);
    assert(isequal(size(L),[1,1]) && abs(L+6)<1e-13, ...
        'For quaternion Hermitian scalar x, A*x+x*A'' must equal 2*Re(A)*x.');
    assert(info.dimension==1,'The scalar Hermitian space has the wrong dimension.');
    [L,info]=qstab.lyapoperator(-2*eye(2));
    assert(isequal(size(L),[6,6]) && norm(L+4*eye(6),2)<1e-13, ...
        'The operator of -2*I must be -4 times the identity on Hermitian space.');
    assert(info.dimension==6,'The 2-by-2 quaternion Hermitian space has dimension six.');
end

function localReentry()
    % Exact Hurwitz condition: alpha>-13/48 and
    % (alpha-1/4)*(alpha-3/4)>0, with the component containing zero.
    [~,R] = qstab.raymargin([13/16,1,1,1],[3,1,1],1);
    localComputed(R);
    localClose(R.positive,1/4,2e-10);
    localClose(R.negative,-13/48,2e-10);
    assert(R.positive<3/4,'The renewed stable region was mistaken for the first exit.');
end

function localSides()
    [P0,P1,u]=localSideData();
    [~,L] = qstab.raymargin(P0,P1,u,'Side','left');
    [~,R] = qstab.raymargin(P0,P1,u,'Side','right');
    localComputed(L); localComputed(R);
    localClose(L.positive,1,2e-10);
    localClose(R.positive,2-sqrt(3),2e-10);
    assert(isinf(L.negative)&&L.negative<0&&isinf(R.negative)&&R.negative<0, ...
        'Both source families are stable on the full negative ray.');
    assert(L.diagnostics.infiniteEigenvalueCount>0 && R.diagnostics.infiniteEigenvalueCount>0, ...
        'This regression must exercise finite contacts together with infinite pencil modes.');
    assert(L.positive<3 && R.positive<2+sqrt(3), ...
        'A later guardian contact was selected instead of the first boundary.');
end

function localInfinity()
    % A pure-imaginary scalar direction changes no real part: L1 is zero.
    A1=localParts(0,1,0,0);
    [~,R] = qstab.raymargin(-1,A1,1,'Type','matrix');
    localComputed(R); localAllRealStable(R);
    [~,R] = qstab.raymargin([1,1],0,1);
    localComputed(R); localAllRealStable(R);

    % An arbitrary strictly upper triangular perturbation preserves diag -1.
    [~,R] = qstab.raymargin(-eye(2),[0,1;0,0],1,'Type','matrix');
    localComputed(R); localAllRealStable(R);

    % A singular direction can also have a genuine finite first contact.
    [~,R] = qstab.raymargin(-eye(2),[1,1;0,0],1,'Type','matrix');
    localComputed(R);
    localClose(R.positive,1,2e-10);
    assert(isinf(R.negative)&&R.negative<0,'The negative endpoint should be infinite.');
end

function localTinyDirection()
    % -1 + alpha*1e-16 first reaches zero at alpha=1e16.
    [~,R] = qstab.raymargin(-1,1e-16,1,'Type','matrix');
    localComputed(R);
    assert(isfinite(R.positive),'A small nonzero direction was silently discarded.');
    localClose(R.positive,1e16,5e-12);
    assert(isinf(R.negative)&&R.negative<0,'The negative scalar endpoint should be infinite.');
end

function localNearReal()
    B=[0,1;-1,0];
    examples={[-1,2-1e-10;0,-1],[-1e8,2e8-1;0,-1e8]};
    % det(A+alpha*B)=alpha^2+c*alpha+a^2 has negative discriminant
    % for both examples, and trace is the fixed negative value -2*a.
    for j=1:numel(examples)
        for tolerance=[1e-9,1e-2]
            [~,R] = qstab.raymargin(examples{j},B,1,'Type','matrix','Tolerance',tolerance);
            assert(~strcmp(R.status,'assumptionFailed'), ...
                'The analytically Hurwitz nominal matrix was rejected.');
            if strcmp(R.positiveStatus,'computed')
                assert(isinf(R.positive)&&R.positive>0, ...
                    'A genuinely complex pencil value became a false positive-side contact.');
            else
                assert(isnan(R.positive),'An unresolved endpoint must remain NaN.');
            end
            if strcmp(R.negativeStatus,'computed')
                assert(isinf(R.negative)&&R.negative<0, ...
                    'A genuinely complex pencil value became a false negative-side contact.');
            else
                assert(isnan(R.negative),'An unresolved endpoint must remain NaN.');
            end
        end
    end
end

function localTangency()
    % p0=s^2+2s+1+j, p1=-(1+j)*(s+1). The first loss is an
    % isolated zero at alpha=1, followed by a renewed loss at alpha=3/2.
    P0=localParts([1,2,1],[0,0,0],[1,0,0],[0,0,0]);
    P1=localParts([-1,-1],[0,0],[-1,-1],[0,0]);
    [~,R] = qstab.raymargin(P0,P1,1);
    if strcmp(R.positiveStatus,'computed')
        % A future improved method may resolve the touching contact.
        localClose(R.positive,1,1e-7);
    else
        assert(isnan(R.positive),'An unresolved touching contact must not become a finite result.');
        assert(any(strcmp(R.status,{'partial','undetermined'})), ...
            'A defective touching contact needs an unresolved/partial status.');
        localClose(R.candidateInterval(2),1,1e-6);
    end
    assert(~(strcmp(R.positiveStatus,'computed') && abs(R.positive-1.5)<1e-4), ...
        'The isolated first loss was skipped in favor of the later simple contact.');
end

function localScaling()
    p=[13/16,1,1,1]; q=[3,1,1];
    for scale=[1e-200,1e200]
        [~,R] = qstab.raymargin(scale*p,scale*q,1);
        localComputed(R);
        localClose(R.positive,1/4,5e-10);
        localClose(R.negative,-13/48,5e-10);
    end
    % A common nonreal quaternion multiplier is legitimate on a real ray.
    P=localParts(2*p,-3*p,4*p,-5*p);
    Q=localParts(2*q,-3*q,4*q,-5*q);
    [~,R] = qstab.raymargin(P,Q,1);
    localComputed(R);
    localClose(R.positive,1/4,5e-9);
    localClose(R.negative,-13/48,5e-9);
end

function localDirection()
    % The supplied direction is u=2, so -1+alpha*u has alpha*=1/2.
    [~,R] = qstab.raymargin(-1,1,2,'Type','matrix');
    localComputed(R);
    localClose(R.positive,1/2,1e-12);
    localClose(R.directionNorm,2,1e-12);
    localClose(R.positiveParameterNorm,1,1e-12);
    assert(isequal(R.directionComponents,[2;0;0;0]), ...
        'The raw direction was silently normalized.');
end

function localAssumptions()
    [~,R] = qstab.raymargin(1,1,1,'Type','matrix');
    assert(strcmp(R.status,'assumptionFailed') && isnan(R.positive), ...
        'A non-Hurwitz nominal member cannot define the requested stable component.');
    [~,R] = qstab.raymargin([1,1],[1,1],1);
    assert(strcmp(R.status,'assumptionFailed'), ...
        'The fixed-degree scalar ray assumption deg(p1)<deg(p0) was not enforced.');
    caught=false;
    try
        qstab.raymargin([1,1],1,0);
    catch
        caught=true;
    end
    assert(caught,'The zero ray direction was accepted.');
end

function [P0,P1,u]=localSideData()
    P0=localParts([1,2,1],[0,0,0],[1,0,0],[0,0,0]);
    P1=localParts([0,0],[0,1],[0,0],[1,0]);
    u=localParts(0,1,0,1);
end

function S=localParts(a,b,c,d)
    S=struct('a',a,'b',b,'c',c,'d',d);
end

function localComputed(R)
    assert(strcmp(R.status,'computed') && strcmp(R.positiveStatus,'computed') && ...
        strcmp(R.negativeStatus,'computed'),'The analytically resolved ray endpoints were not computed.');
end

function localAllRealStable(R)
    assert(isinf(R.negative)&&R.negative<0&&isinf(R.positive)&&R.positive>0, ...
        'A family stable on the whole real line acquired a false finite endpoint.');
end

function localClose(actual,expected,tolerance)
    assert(isfinite(actual) && abs(actual-expected)<=tolerance*max(1,abs(expected)), ...
        'The computed value disagrees with its analytic reference.');
end
