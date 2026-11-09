function tests = qstab_cases_balls_native()
tests={ ...
 'Native quaternion coefficient-ball distances agree with component arrays',@nativeValues; ...
 'Native quaternion CT and Schur balls preserve full class directions',@nativeMargins; ...
 'Native quaternion critical members leave fixed coefficients unchanged',@nativeWitness; ...
 'Native delay coefficient balls preserve grade and reversed weights',@nativeDelay};
end
function nativeValues()
a=[2 .5 1];b=[.1 .2 0];c=[.3 -.1 0];d=[-.2 .1 0];P=quaternion(a,b,c,d);
Q=struct('a',a,'b',b,'c',c,'d',d);
[r,I]=qstab.ballvalue(P,[1 .5 0],[0 1i 2+.1i],'Witness',true);s=qstab.ballvalue(Q,[1 .5 0],[0 1i 2+.1i]);
assert(all(I.valid)&&max(abs(r-s))<1e-12);
end
function nativeMargins()
P=quaternion([2 1],[0 0],[3 0],[0 0]);[~,R] = qstab.ballradius(P,[1 0], ...
    'AbsTol',2e-3,'RelTol',2e-3,'MaxCells',2047,'MaxTime',Inf);
assert(strcmp(R.status,'converged')&&R.lower<=2&&R.upper>=2-1e-12);
P=quaternion([-.3 1],[0 0],[-.4 0],[0 0]);[~,R] = qstab.ballradius(P,[1 0], ...
    'Variable','z','AbsTol',2e-3,'RelTol',2e-3,'MaxTime',Inf);
assert(strcmp(R.status,'converged')&&R.lower<=.5&&R.upper>=.5-1e-12);
end
function nativeWitness()
P=quaternion([2 3 1],[.1 .2 0],[.3 -.1 0],[0 0 0]);C=qstab.polycomponents(P);
[~,I]=qstab.ballvalue(P,[1 0 0],1i,'Witness',true);W=I.witnesses{1};
assert(I.valid&&W.valid&&isequal(C(:,2:3),W.coefficientComponents(:,2:3)));
N=quaternion(reshape(W.coefficientComponents(1,:),1,1,[]), ...
    reshape(W.coefficientComponents(2,:),1,1,[]),reshape(W.coefficientComponents(3,:),1,1,[]), ...
    reshape(W.coefficientComponents(4,:),1,1,[]));
[lambda,info]=qstab.roots(N);assert(info.isValid&&min(abs(lambda-1i))<1e-5);
end
function nativeDelay()
P=quaternion([1 -.5 0],[0 0 0],[0 0 0],[0 0 0]);
[~,R] = qstab.ballradius(P,[0 1 0],'Variable','zinv','AbsTol',2e-3,'RelTol',2e-3,'MaxTime',Inf);
assert(strcmp(R.status,'converged')&&R.grade==2&&R.lower<=.5&&R.upper>=.5-1e-12);
assert(isequal(R.witness.weights,[0 1 0])&&all(R.witness.coefficientComponents(:,3)==0));
end
