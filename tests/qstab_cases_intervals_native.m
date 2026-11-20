function cases=qstab_cases_intervals_native()
%QSTAB_CASES_INTERVALS_NATIVE Real MATLAB quaternion input, no emulation.
cases={'Native quaternion inputs reproduce all current interval examples',@allNative; ...
 'Native single-component frequency witness is an actual non-i right zero',@witness};
end
function allNative()
R=run_qstab_intervals('UseQuaternion',true,'Save',false,'Verbose',false);
assert(R.ok&&R.nativeQuaternionInputs&&~R.isOctave&&numel(R.cases)==24);
end
function witness()
P=quaternion([1 1 1],[2 0 0],[0 2 0],[0 2 0]);W=zeros(4,3);W(1,1:2)=.1;
u=fzero(@(x)x^3+8*x-4,[4/9,.5]);
[r,I]=qstab.boxfreqmargin(P,W,sqrt(u),'Witness',true);
assert(I.valid&&abs(r-10*(1-u))<1e-9);
Z=I.witnesses{1};assert(Z.valid&&abs(Z.boundaryComponents(3))>.1&&abs(Z.boundaryComponents(4))>.1);
end
