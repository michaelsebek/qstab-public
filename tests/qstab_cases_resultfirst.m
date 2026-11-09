function cases=qstab_cases_resultfirst()
%QSTAB_CASES_RESULTFIRST Exercise the public entry points as an application.
cases={
 'A single decision output is a double; complete diagnostics remain available',@nominal
 'Strict boundary, unresolved band and invalid grade are distinguishable',@boundary
 'All five decision functions return the same scalar convention',@allTests
 'Algebraic norm values and Bernstein regions can be used as primary outputs',@sets
 'Stored certificates and actual witness coefficients can be reused',@certificates
 'Signed ray endpoints and global brackets keep their different meanings',@margins
 'Spectral-set membership preserves the query layout and undecided points',@membership
};
end
function nominal()
p=[2 3 1]; t=qstab.isstable(p); [u,I]=qstab.isstable(p);
assert(isa(t,'double')&&isscalar(t)&&t==1&&u==t&&I.isStable);
assert(isfield(I,'diagnostics')&&isfield(I,'method')&&~I.spectrumComputed);
[t,I]=qstab.isstable([-1 1]); assert(t==0&&strcmp(I.status,'notStable'));
[t,I]=qstab.isstable(0); assert(isnan(t)&&strcmp(I.status,'assumptionFailed'));
end
function boundary()
[t,I]=qstab.isstable([0 1]); assert(t==0&&strcmp(I.status,'notStable')&&~I.spectrumComputed);
[t,I]=qstab.isstable([1i 1]); assert(t==0&&strcmp(I.status,'notStable'));
[t,I]=qstab.isstable([-1+1e-12 1],'Variable','z','Tolerance',1e-8);
assert(isnan(t)&&strcmp(I.status,'undetermined'));
[t,I]=qstab.isstable([0 1],'Variable','zinv'); assert(isnan(t)&&strcmp(I.status,'assumptionFailed'));
end
function allTests()
p=[2 1]; V={p,[3 1]}; W=zeros(4,2); W(1,1)=1;
for entry={@()qstab.isstable(p),@()qstab.robuststable(V),@()qstab.zeroexclude(V), ...
           @()qstab.ballstable(p,[1 0],0),@()qstab.boxstable(p,W,0)}
 t=entry{1}();assert(isa(t,'double')&&isscalar(t)&&t==1);
end
[t,I]=qstab.robuststable({[-1 1]});assert(t==0&&strcmp(I.status,'notStable')&&~isempty(I.witness));
end
function sets()
V={[1 1],[3 1]}; [v,I]=qstab.valueset(V,1i,'Weights',[.5 .5]);assert(abs(v-(2+1i)^2)<1e-12&&isequal(v,I.points));
[E,J]=qstab.valueenclosure(V,1i);assert(isequaln(E,J.regions)&&~isempty(E));
[d,D]=qstab.leadingbound(V);assert(d==1&&D.isSeparated);
end
function certificates()
V={[1 1],[2 1]};[t,R]=qstab.zeroexclude(V,'StoreCertificates',true);
[rows,D]=qstab.certificateData(R);assert(t==1&&~isempty(rows)&&isequaln(rows,D.rows));
[regions,E]=qstab.certificateEnclosure(V,R);assert(~isempty(regions)&&isequaln(regions,E.regions));
[C,W]=qstab.witnessData(V,R);assert(isempty(C)&&~W.available);
V={[-1 1]};[t,R]=qstab.robuststable(V);[C,W]=qstab.witnessData(V,R);
assert(t==0&&W.isAdmissible&&isequal(C,[-1 1;zeros(3,2)]));
end
function margins()
p=[2 1]; [ab,R]=qstab.raymargin(p,1,1);assert(isequaln(ab,R.interval)&&abs(ab(1)+2)<1e-10&&isinf(ab(2)));
[b,I]=qstab.stabilityradius(p,1);assert(isequaln(b,[I.lower I.upper])&&b(1)<=2+1e-10&&b(2)>=2-1e-10);
[b,I]=qstab.ballradius(p,[1 0]);assert(isequaln(b,[I.lower I.upper])&&b(1)<=2+1e-10&&b(2)>=2-1e-10);
W=zeros(4,2);W(1,1)=1;[b,I]=qstab.boxradius(p,W,'MaxIterations',2);assert(isequaln(b,[I.lower I.upper])&&b(1)<=2+1e-10&&b(2)>=2-1e-10);
end
function membership()
[x,I]=qstab.ballset([2 1],[1 0],[-2 -1 0],1);
assert(isequaln(x,[1 NaN 0])&&I.inside(1)&&I.unresolved(2)&&I.outside(3));
end
