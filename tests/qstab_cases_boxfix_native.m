function tests=qstab_cases_boxfix_native()
tests={ ...
 'Native quaternion unstable CT box returns an original-box spectral witness',@continuous; ...
 'Native delay box corner keeps the zero forward mode and fixed original coefficients',@delay ...
};
end
function continuous()
P=quaternion([2 1],[0 0],[3 0],[0 0]);W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxstable(P,W,3,'MaxTime',Inf);
assert(strcmp(R.status,'notStable')&&R.witness.valid&&R.witness.insideRequestedBox);
assert(R.witness.weightedNorm<=3&&R.witness.coefficientComponents(3,1)==3);
assert(R.witness.spectralCheck.isConfirmedCounterexample);
end
function delay()
P=quaternion([1 -.5 0],[0 0 0],[0 0 0],[0 0 0]);W=zeros(4,3);W(1,2)=1;
[~,R] = qstab.boxstable(P,W,1,'Variable','zinv','MaxTime',Inf);
assert(strcmp(R.status,'notStable')&&R.grade==2&&R.witness.insideRequestedBox);
assert(R.witness.weightedNorm<=1&&all(R.witness.coefficientComponents(:,3)==0));
assert(R.witness.coefficientComponents(1,1)==1);
assert(any(R.witness.spectralCheck.testRepresentatives==0));
assert(any(isinf(R.witness.spectralCheck.inputRepresentatives)));
end
