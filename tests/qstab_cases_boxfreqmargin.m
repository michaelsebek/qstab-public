function tests=qstab_cases_boxfreqmargin()
%QSTAB_CASES_BOXFREQMARGIN Two-point class margins for one Hamilton component.
tests={ ...
 'One-component slice returns 3-sqrt(5) and an actual non-i right root',@slice; ...
 'Two independent regulator gains reproduce the algebraic critical frequency',@regulator; ...
 'Even and odd widths remain independent and the maximum of their costs is used',@parities; ...
 'Original coefficient lifting in all four components agrees with the independent point gauge',@liftAllComponents; ...
 'Structurally absent parity gives infinity without tolerance projection',@frozenParity; ...
 'Zero widths distinguish an existing class zero from unattainable contact',@zeroWidths; ...
 'A zero nominal polynomial gives zero without a global stability assertion',@zeroNominal; ...
 'Grade padding and row or column frequency shape are preserved',@gradeAndShape; ...
 'Extreme common coefficient and width scaling preserves the class margin',@commonScaling; ...
 'Unrepresentable positive radii and underflowed powers never become zero or infinity claims',@rangeGuards; ...
 'Near-cancelled frozen parity remains undetermined and does not become a contact',@nearFrozen; ...
 'Invalid frequencies, widths, variables and multiple uncertain components fail',@invalid ...
};
end

function slice()
C=[1 2 1;0 3 0;3 0 0;0 0 0];W=zeros(size(C));W(2,2)=1;
[r,I]=qstab.boxfreqmargin(parts(C),W,1,'Witness',true);
assert(I.valid&&abs(r-(3-sqrt(5)))<1e-12&&I.componentIndex==2);
assert(~I.globalMinimumClaim&&~I.validated&&isequal(size(I.branchRadii),[2 1]));
assert(all(abs(sort(I.branchRadii)-[3-sqrt(5);3+sqrt(5)])<1e-12));
V=I.witnesses{1};checkWitness(V,C,W,r);
assert(abs(V.coefficientComponents(2,2)-sqrt(5))<1e-12);
assert(abs(V.boundaryComponents(3)+2/3)<1e-12);
assert(abs(V.boundaryComponents(4)-sqrt(5)/3)<1e-12);
assert(abs(V.boundaryComponents(2))<1e-12);
% Query without retention still performs reconstruction, but keeps no payload.
[r2,J]=qstab.boxfreqmargin(parts(C),W,1);
assert(abs(r2-r)<eps&&isempty(J.witnesses{1}));
end

function regulator()
u=roots([1 0 8 -4]);u=real(u(abs(imag(u))<1e-12));
C=[1 1 1;2 0 0;0 2 0;0 2 0];W=zeros(size(C));W(1,1:2)=.1;
[r,I]=qstab.boxfreqmargin(parts(C),W,[0 sqrt(u) 1],'Witness',true);
assert(all(I.valid)&&isinf(r(1))&&abs(r(2)-10*(1-u))<2e-12&&abs(r(3)-20)<2e-12);
V=I.witnesses{2};checkWitness(V,C,W,r(2));
assert(max(abs(V.coefficientComponents(1,1:2)-u))<1e-12);
assert(norm(V.boundaryComponents-[0;-u*u/2;u;-u])<2e-12);
assert(abs(r(2)^3-30*r(2)^2+1100*r(2)-5000)<2e-10);
end

function parities()
W=zeros(4,2);W(1,:)=[2 .5];C=[2 1;zeros(3,2)];
[r,I]=qstab.boxfreqmargin([2 1],W,1,'Witness',true);
assert(I.valid&&abs(r-2)<1e-12);
assert(max(abs(I.frequencyData{1}.branches(1).parityRadii-[1;2]))<1e-12);
checkWitness(I.witnesses{1},C,W,r);
end


function liftAllComponents()
C=[2 3 4 2 1;1 -2 .3 .4 0;.2 .5 -.8 0 .1;1 .2 0 -.3 .2];
for nu=1:4
    W=zeros(size(C));W(nu,:)=[.7 .3 .6 .5 .2];
    [r,I]=qstab.boxfreqmargin(parts(C),W,[.4 1.3],'Witness',true);
    assert(all(I.valid));
    for k=1:2
        V=I.witnesses{k};checkWitness(V,C,W,r(k));
        point=struct('a',V.boundaryComponents(1),'b',V.boundaryComponents(2), ...
            'c',V.boundaryComponents(3),'d',V.boundaryComponents(4));
        [pointRadius,J]=qstab.boxvalue(parts(C),W,point);
        assert(J.valid&&abs(pointRadius-r(k))<1e-9*max(r(k),1));
    end
end
end

function frozenParity()
W=zeros(4,2);W(1,1)=.1;
[r,I]=qstab.boxfreqmargin([1 1],W,[0 1],'Tolerance',.01,'Witness',true);
assert(all(I.valid)&&abs(r(1)-10)<1e-12&&isinf(r(2)));
C=[1+1e-8 1 1;zeros(3,3)];W=zeros(size(C));W(1,2)=1;
[r,I]=qstab.boxfreqmargin(parts(C),W,1,'Tolerance',.1);
assert(I.valid&&isinf(r));
end

function zeroWidths()
C=[0 1;0 0;-1 0;0 0];W=zeros(size(C));
[r,I]=qstab.boxfreqmargin(parts(C),W,[0 1 2],'Witness',true);
assert(all(I.valid)&&isinf(r(1))&&r(2)==0&&isinf(r(3))&&I.componentIndex==0);
assert(isequal(I.witnesses{2}.coefficientComponents,C));
assert(norm(I.witnesses{2}.boundaryComponents-[0;0;1;0])<1e-12);
end

function zeroNominal()
W=zeros(4,3);W(4,1)=1;
[r,I]=qstab.boxfreqmargin([0 0 0],W,[0 1 100],'Witness',true);
assert(all(r==0)&&all(I.valid)&&~I.globalMinimumClaim);
assert(all(cellfun(@(V)V.weightedNorm==0,I.witnesses)));
end

function gradeAndShape()
W=zeros(4,4);W(1,1)=1;
[r,I]=qstab.boxfreqmargin([1 1],W,[0;2],'Grade',3,'Witness',true);
assert(isequal(size(r),[2 1])&&I.grade==3&&all(I.valid)&&abs(r(1)-1)<1e-12&&isinf(r(2)));
assert(isequal(size(I.witnesses{1}.coefficientComponents),[4 4]));
end

function commonScaling()
C=[1 2 1;0 3 0;3 0 0;0 0 0];W=zeros(size(C));W(2,2)=1;
for scale=[1e-200 1 1e200]
    [r,I]=qstab.boxfreqmargin(parts(scale*C),scale*W,1,'Witness',true);
    assert(I.valid&&abs(r-(3-sqrt(5)))<2e-12);
    checkWitness(I.witnesses{1},scale*C,scale*W,r);
end
end

function rangeGuards()
W=zeros(4,2);W(1,1)=1e-300;
[r,I]=qstab.boxfreqmargin(1e300*[1 1],W,0,'Witness',true);
assert(isnan(r)&&~I.valid&&isempty(I.witnesses{1}));
W(1,1)=1e300;
[r,I]=qstab.boxfreqmargin(1e-300*[1 1],W,0,'Witness',true);
assert(isnan(r)&&~I.valid&&isempty(I.witnesses{1}));
W=zeros(4,3);W(1,1)=1;
[r,I]=qstab.boxfreqmargin([1 0 1],W,1e-200);
assert(isnan(r)&&~I.valid);
end

function nearFrozen()
W=zeros(4,3);W(1,2)=1;
[r,I]=qstab.boxfreqmargin([1+eps 1 1],W,1,'Tolerance',.01);
assert(isnan(r)&&~I.valid);
[r,I]=qstab.boxfreqmargin([1 1 1],W,1,'Witness',true);
assert(I.valid&&abs(r-1)<1e-12);
end

function invalid()
W=zeros(4,2);W(1,1)=1;
fail(@()qstab.boxfreqmargin([1 1],W,-1),'qstab:boxfreqmargin:Frequency');
fail(@()qstab.boxfreqmargin([1 1],W,Inf),'qstab:boxfreqmargin:Frequency');
fail(@()qstab.boxfreqmargin([1 1],W,1i),'qstab:boxfreqmargin:Frequency');
fail(@()qstab.boxfreqmargin([1 1],W,[]),'qstab:boxfreqmargin:Frequency');
fail(@()qstab.boxfreqmargin([1 1],W,1,'Variable','z'),'qstab:boxfreqmargin:Variable');
fail(@()qstab.boxfreqmargin([1 1],parts(W),1),'qstab:boxfreqmargin:Weights');
fail(@()qstab.boxfreqmargin([1 1],-W,1),'qstab:box:Weights');
fail(@()qstab.boxfreqmargin([1 1],zeros(4,3),1),'qstab:box:Weights');
W(2,2)=1;
fail(@()qstab.boxfreqmargin([1 1],W,1),'qstab:boxfreqmargin:Components');
end

function checkWitness(V,C,W,r)
assert(V.valid&&V.fixedComponentsPreserved&&~V.validated);
assert(V.evaluationResidual<1e-10&&abs(V.weightedNorm-r)<1e-10*max(1,r));
D=V.coefficientComponents-C;active=W>0;
assert(all(D(~active)==0));
assert(all(abs(D(active))./W(active)<=r*(1+1e-10)));
assert(isequal(D,V.perturbationComponents));
end
function fail(f,id)
ok=false;
try,f();catch problem,ok=strcmp(problem.identifier,id);end
assert(ok);
end
function P=parts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
