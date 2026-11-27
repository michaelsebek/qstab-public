function L=exactOneParameterLifting(C,N,g,s,h,lo,hi)
%EXACTONEPARAMETERLIFTING Contact verification in a selected real field.
% The real algebraic gamma is represented by a square-free modulus and an
% isolating interval. Arithmetic is reduced modulo that modulus. Whenever
% a zero divisor occurs, exact Sturm counts retain the factor containing
% the selected root. No approximate root is used in a field decision.
L=struct('status','notComputed','exact',false,'numericalWitness',false, ...
    'method','original-polynomial-over-selected-real-algebraic-field', ...
    'gammaExactBounds',{{char(lo),char(hi)}},'fieldModulus','','fieldVariable',char(g), ...
    'frequencyVariable','qstabExactX', ...
    'zeroContact',false,'positiveFrequencyCount',0,'frequencySquared','', ...
    'frequencyPolynomialCoefficients',{{}},'remainderA',{{}},'remainderB',{{}}, ...
    'rightEvaluationZeroComponents',{{}},'leftEvaluationZeroComponents',{{}}, ...
    'rightResidualExactlyZero',false,'leftResidualExactlyZero',false, ...
    'sphereOfEvaluationZeros',false,'normIdentitiesVerified',false,'message','');
rationalRoot=isequal(lo,hi);
if rationalRoot,m=g-lo;else,[~,m]=qstab.internal.exactOneParameter('sturm',h,g);end
gammaLo=lo;gammaHi=hi;
if rationalRoot,C=subs(C,g,lo);N=subs(N,g,lo);end
for k=1:numel(C),C(k)=red(C(k));end
if allZeroField(C(:,1))
    L.zeroContact=true;L.status='verified';L.exact=true;
    L.rightEvaluationZeroComponents={'0';'0';'0';'0'};
    L.leftEvaluationZeroComponents={'0';'0';'0';'0'};
    L.rightResidualExactlyZero=true;L.leftResidualExactlyZero=true;
    L.message='The original quaternion constant coefficient vanishes exactly at the first contact.';
    finish();return
end
a=fliplr(coeffs(expand(N),s,'All'));e={};o={};
for k=0:numel(a)-1
    if mod(k,2)==0,e{end+1}=red((-1)^(k/2)*a(k+1)); %#ok<AGROW>
    else,o{end+1}=red((-1)^((k-1)/2)*a(k+1));end %#ok<AGROW>
end
e=trim(e);o=trim(o);D=pgcd(e,o);
if numel(D)<2
    L.status='inconsistent';L.message='No imaginary-axis class could be lifted from the exact first-root candidate.';
    finish();return
end
derivative=cell(1,numel(D)-1);
for k=2:numel(D),derivative{k-1}=red((k-1)*D{k});end
D=pdivide(D,pgcd(D,derivative));
x=sym('qstabExactX');L.frequencyPolynomialCoefficients=textCells(D);
if numel(D)==2
    xx=red(-D{1}/D{2});
    if signField(xx)<=0
        L.status='inconsistent';L.message='The common norm factor has no positive squared frequency.';finish();return
    end
    L.positiveFrequencyCount=1;L.frequencySquared=char(xx);
    [A,B]=remainder(xx);b2=red(sum(B.^2));
    L.remainderA=qstab.internal.exactOneParameter('strings',A);
    L.remainderB=qstab.internal.exactOneParameter('strings',B);
    L.normIdentitiesVerified=zero(red(sum(A.^2)-xx*b2))&&zero(red(sum(A.*B)));
    if ~L.normIdentitiesVerified,error('qstab:exact:Internal','Original-polynomial norm identities failed during exact lifting.');end
    if zero(b2)
        if ~allZeroField(A)||~allZeroField(B),error('qstab:exact:Internal','A zero quaternion norm had nonzero real components.');end
        L.sphereOfEvaluationZeros=true;L.status='verified';L.exact=true;
        L.rightResidualExactlyZero=true;L.leftResidualExactlyZero=true;
        L.message='The central remainder is zero: the whole imaginary sphere consists of evaluation zeros.';
    else
        bc=[B(1);-B(2:4)];
        lr=-qstab.internal.exactOneParameter('mul',bc,A)/b2;
        ll=-qstab.internal.exactOneParameter('mul',A,bc)/b2;
        for k=1:4,lr(k)=red(lr(k));ll(k)=red(ll(k));end
        L.rightEvaluationZeroComponents=qstab.internal.exactOneParameter('strings',lr);
        L.leftEvaluationZeroComponents=qstab.internal.exactOneParameter('strings',ll);
        vr=C(:,end);vl=vr;
        for k=size(C,2)-1:-1:1
            vr=qstab.internal.exactOneParameter('mul',vr,lr)+C(:,k);
            vl=qstab.internal.exactOneParameter('mul',ll,vl)+C(:,k);
            for j=1:4,vr(j)=red(vr(j));vl(j)=red(vl(j));end
        end
        L.rightResidualExactlyZero=allZeroField(vr);L.leftResidualExactlyZero=allZeroField(vl);
        if ~L.rightResidualExactlyZero||~L.leftResidualExactlyZero
            error('qstab:exact:Internal','An exact reconstructed evaluation zero failed on the original polynomial.');
        end
        L.status='verified';L.exact=true;
        L.message='Both evaluation zeros and their original-polynomial residuals were verified in the selected real algebraic field.';
    end
else
    % A general algebraic contact remains a triangular exact representation:
    % m(gamma)=0, D(x,gamma)=0, x>0. Sturm over the ordered real gamma-field
    % proves that such a class exists; central remainder identities lift it.
    DS=psturm(D);L.positiveFrequencyCount=pvariation(DS,sym(0))-pvariation(DS,Inf);
    if L.positiveFrequencyCount<1
        L.status='inconsistent';L.message='The exact common norm factor has no positive real squared frequency.';finish();return
    end
    [A,B]=remainder(x);L.remainderA=qstab.internal.exactOneParameter('strings',A);
    L.remainderB=qstab.internal.exactOneParameter('strings',B);
    u=polyCoefficients(expand(sum(A.^2)-x*sum(B.^2)),x);
    v=polyCoefficients(expand(sum(A.*B)),x);
    L.normIdentitiesVerified=polyZero(prem(u,D))&&polyZero(prem(v,D));
    if ~L.normIdentitiesVerified,error('qstab:exact:Internal','The general algebraic lifting identities failed.');end
    L.status='verifiedImplicit';L.exact=true;
    L.frequencySquared='positive real roots of frequencyPolynomialCoefficients in the selected gamma field';
    L.message=['An exact positive-frequency class is proved by ordered-field Sturm counting. ', ...
        'Its remainder satisfies |A|^2=x|B|^2 and Re(A*conj(B))=0. ', ...
        'For B~=0 use -inv(B)*A (right evaluation) or -A*inv(B) (left); ', ...
        'for B=0 the whole class consists of zeros. No numerical witness is asserted.'];
end
finish();

    function finish()
        L.fieldModulus=char(m);L.gammaExactBounds={char(gammaLo),char(gammaHi)};
    end
    function [A,B]=remainder(xx)
        A=sym(zeros(4,1));B=A;
        for j=0:size(C,2)-1
            if mod(j,2)==0,A=A+C(:,j+1)*(-xx)^(j/2);
            else,B=B+C(:,j+1)*(-xx)^((j-1)/2);end
        end
        free=symvar(xx);dependsOnX=false;
        for jj=1:numel(free),dependsOnX=dependsOnX||isequal(free(jj),x);end
        if ~dependsOnX
            for j=1:4,A(j)=red(A(j));B(j)=red(B(j));end
        else
            A=expand(A);B=expand(B);
        end
    end
    function tf=zero(a)
        a=red(a);
        if qstab.internal.exactOneParameter('iszero',a),tf=true;return;end
        if rationalRoot||isempty(symvar(a)),tf=false;return;end
        [nu,~]=numden(a);c=qstab.internal.exactOneParameter('gcd',m,nu,g);
        if numel(coeffs(c,g,'All'))<2,tf=false;return;end
        if contains(c),m=c;tf=true;else,[m,r]=quorem(m,c,g);m=expand(m);requireZero(r);tf=false;end
    end
    function a=red(a)
        if isempty(symvar(a)),a=simplify(a);return;end
        if rationalRoot,a=simplify(subs(a,g,lo));return;end
        [nu,de]=numden(simplify(a));
        if isempty(symvar(de))
            [~,a]=quorem(nu,m,g);a=expand(a/de);return
        end
        [~,nu]=quorem(nu,m,g);[~,de]=quorem(de,m,g);
        de=expand(de);nu=expand(nu);
        c=qstab.internal.exactOneParameter('gcd',m,de,g);
        if numel(coeffs(c,g,'All'))>1
            if contains(c),error('qstab:exact:Internal','A selected algebraic-field denominator vanishes.');end
            [m,r]=quorem(m,c,g);m=expand(m);requireZero(r);
            [~,nu]=quorem(nu,m,g);[~,de]=quorem(de,m,g);
        end
        inv=inverse(de);[~,a]=quorem(expand(nu*inv),m,g);a=expand(a);
    end
    function v=inverse(a)
        r0=m;r1=a;t0=sym(0);t1=sym(1);
        while ~qstab.internal.exactOneParameter('iszero',r1)
            [q,r]=quorem(r0,r1,g);r0=r1;r1=expand(r);
            t=t0-q*t1;t0=t1;t1=expand(t);
        end
        if numel(coeffs(r0,g,'All'))~=1,error('qstab:exact:Internal','Algebraic-field inversion is not coprime.');end
        [~,v]=quorem(t0/r0,m,g);v=expand(v);
    end
    function tf=contains(p)
        [S,~]=qstab.internal.exactOneParameter('sturm',p,g);
        tf=qstab.internal.exactOneParameter('count',S,g,gammaLo,gammaHi)>0;
    end
    function v=signField(a)
        if zero(a),v=0;return;end
        a=red(a);
        if rationalRoot||isempty(symvar(a)),v=qstab.internal.exactOneParameter('sign',a);return;end
        [S,~]=qstab.internal.exactOneParameter('sturm',a,g);
        for z=1:512
            va=qstab.internal.exactOneParameter('sign',subs(a,g,gammaLo));
            if va~=0&&qstab.internal.exactOneParameter('count',S,g,gammaLo,gammaHi)==0,v=va;return;end
            mid=(gammaLo+gammaHi)/2;
            [MS,~]=qstab.internal.exactOneParameter('sturm',m,g);
            if qstab.internal.exactOneParameter('count',MS,g,gammaLo,mid)>0,gammaHi=mid;else,gammaLo=mid;end
        end
        error('qstab:exact:Undecidable','The ordered algebraic-field sign refinement budget was exhausted.');
    end
    function tf=allZeroField(a)
        tf=true;for z=1:numel(a),if ~zero(a(z)),tf=false;return;end,end
    end
    function p=trim(p)
        if isempty(p),p={sym(0)};return;end
        while numel(p)>1&&zero(p{end}),p(end)=[];end
        for z=1:numel(p),p{z}=red(p{z});end
    end
    function tf=polyZero(p),p=trim(p);tf=numel(p)==1&&zero(p{1});end
    function r=prem(a,b)
        a=trim(a);b=trim(b);
        if polyZero(b),error('qstab:exact:Internal','Polynomial division by zero in algebraic lifting.');end
        while numel(a)>=numel(b)&&~polyZero(a)
            k=numel(a)-numel(b);q=red(a{end}/b{end});
            for z=1:numel(b),a{z+k}=red(a{z+k}-q*b{z});end
            a=trim(a);
        end
        r=a;
    end
    function a=pgcd(a,b)
        a=trim(a);b=trim(b);
        while ~polyZero(b),r=prem(a,b);a=b;b=r;end
        lead=a{end};for z=1:numel(a),a{z}=red(a{z}/lead);end
    end
    function S=psturm(p)
        d=cell(1,max(1,numel(p)-1));
        if numel(p)==1,d={sym(0)};else,for z=2:numel(p),d{z-1}=red((z-1)*p{z});end,end
        c=pgcd(p,d);f=pdivide(p,c);
        if numel(f)==1,S={f};return;end
        d=cell(1,numel(f)-1);for z=2:numel(f),d{z-1}=red((z-1)*f{z});end
        S={f,d};
        while true
            r=prem(S{end-1},S{end});if polyZero(r),break;end
            for z=1:numel(r),r{z}=-r{z};end
            S{end+1}=r; %#ok<AGROW>
        end
    end
    function q=pdivide(a,b)
        a=trim(a);b=trim(b);q=repmat({sym(0)},1,numel(a)-numel(b)+1);
        while numel(a)>=numel(b)&&~polyZero(a)
            k=numel(a)-numel(b);v=red(a{end}/b{end});q{k+1}=v;
            for z=1:numel(b),a{z+k}=red(a{z+k}-v*b{z});end
            a=trim(a);
        end
        if ~polyZero(a),error('qstab:exact:Internal','An exact algebraic polynomial division left a remainder.');end
        q=trim(q);
    end
    function v=pvariation(S,a)
        signs=zeros(1,numel(S));
        for z=1:numel(S)
            p=S{z};
            if isnumeric(a)&&isinf(a),value=p{end};else
                value=p{end};for y=numel(p)-1:-1:1,value=red(value*a+p{y});end
            end
            signs(z)=signField(value);
        end
        signs(signs==0)=[];v=sum(signs(1:end-1).*signs(2:end)<0);
    end
    function p=polyCoefficients(a,v)
        c=fliplr(coeffs(expand(a),v,'All'));p=cell(1,numel(c));
        for z=1:numel(c),p{z}=red(c(z));end
        p=trim(p);
    end
    function c=textCells(p),c=cellfun(@char,p,'UniformOutput',false);end
    function requireZero(a)
        if ~qstab.internal.exactOneParameter('iszero',a),error('qstab:exact:Internal','An exact factor division failed.');end
    end
end
