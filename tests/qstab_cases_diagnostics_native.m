function tests = qstab_cases_diagnostics_native()
% Native tests append to the existing native group; no private qfiles used.
tests={ ...
 'Native quaternion Bernstein enclosures agree with component representations',@nativePoly; ...
 'Native quaternion interval enclosures preserve component bounds',@nativeInterval; ...
 'Native quaternion witness reconstruction retains declared delay grade',@nativeWitness};
end
function nativePoly()
a=quaternion([1 1],[0 0],[-2 0],[0 0]);b=quaternion([1 1],[0 0],[2 0],[0 0]);
[~,E] = qstab.valueenclosure({a,b},0);assert(min(real(E.regions{1}.controls))==-3);
[~,R] = qstab.zeroexclude({a,b},'StoreCertificates',true,'MaxTime',Inf);
[~,C] = qstab.certificateEnclosure({a,b},R);assert(C.geometryResidual<1e-7);
end
function nativeInterval()
L=quaternion([1 2],zeros(1,2),zeros(1,2),zeros(1,2));
U=quaternion([3 4],zeros(1,2),zeros(1,2),zeros(1,2));
[~,E] = qstab.valueenclosure({L,U},1i,'Uncertainty','interval');
assert(isequal(E.regions{1}.valueBounds,[-15 5;4 24]));
end
function nativeWitness()
a=quaternion([1 -1.2 0],zeros(1,3),zeros(1,3),zeros(1,3));
b=quaternion([1 -1.4 0],zeros(1,3),zeros(1,3),zeros(1,3));
[~,R] = qstab.zeroexclude({a,b},'Variable','zinv','MaxTime',Inf);[~,W] = qstab.witnessData({a,b},R);
assert(W.isConfirmedCounterexample&&W.grade==2&&sum(isinf(abs(W.inputRepresentatives)))==1);
end
