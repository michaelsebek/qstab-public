function tests=qstab_cases_boxes_native()
tests={ ...
 'Native quaternion actual-point box gauges agree with component data',@point; ...
 'Native quaternion point arrays retain non-i directions and shape',@arrays; ...
 'Native zero-point box witnesses leave fixed components unchanged',@zero; ...
 'Native component boxes retain CT global radius bounds',@ct; ...
 'Native forward and delay boxes retain grade and the same physical margin',@delay; ...
 'Native qB02 edge witness retains its noncommuting critical direction',@edge ...
};
end
function point()
P=quaternion([2 1],[.3 0],[-.4 0],[.2 0]);W=ones(4,2);W(:,end)=0;
[r,I]=qstab.boxvalue(P,W,quaternion(0,1,0,0),'Witness',true);
assert(I.valid&&isfinite(r)&&I.witnesses{1}.evaluationResidual<1e-8);
end
function arrays()
P=quaternion([2 1],[0 0],[3 0],[0 0]);W=ones(4,2);W(:,end)=0;
q=quaternion(zeros(2),[1 0;0 0],[0 1;0 0],[0 0;1 0]);[r,I]=qstab.boxvalue(P,W,q);
assert(isequal(size(r),[2 2])&&all(I.valid(:)));
end
function zero()
P=quaternion([2 1],[3 .2],[4 -.3],[5 .1]);W=ones(4,2);
[~,I]=qstab.boxvalue(P,W,quaternion(0,0,0,0),'Witness',true);Q=I.witnesses{1};
assert(I.valid&&all(Q.coefficientComponents(:,1)==0)&&Q.evaluationResidual==0);
end
function ct()
P=quaternion([2 1],[0 0],[0 0],[0 0]);W=zeros(4,2);W(1,1)=1;
[~,R] = qstab.boxradius(P,W,'RelTol',.02,'MaxCells',2048,'MaxTime',Inf);
assert(R.lower<=2&&R.upper>=2-1e-9&&R.lower>0);
end
function delay()
P=quaternion([1 -.5 0],[0 0 0],[0 0 0],[0 0 0]);W=zeros(4,3);W(1,2)=1;
[~,R] = qstab.boxradius(P,W,'Variable','zinv','RelTol',.03,'MaxCells',2048,'MaxTime',Inf);
assert(R.lower<=.5&&R.upper>=.5-1e-9&&R.grade==2&&R.witness.valid);
end
function edge()
P=quaternion([1 1 1],[2 0 0],[0 2.5 0],[0 0 0]);W=zeros(4,3);W(3:4,2)=1;
q=quaternion(0,-.5,0,-sqrt(3)/2);[r,I]=qstab.boxvalue(P,W,q,'Witness',true);
assert(I.valid&&abs(r-(2.5-sqrt(3)))<1e-8&&I.witnesses{1}.valid);
end
