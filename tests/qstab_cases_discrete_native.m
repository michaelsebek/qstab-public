function tests = qstab_cases_discrete_native()
%QSTAB_CASES_DISCRETE_NATIVE Genuine MathWorks quaternion input checks.
tests={ ...
 'Native quaternion Schur rays agree with numeric coefficient rays',@ray; ...
 'Native noncommuting unit-class witnesses retain the multiplication side',@quotients; ...
 'Native quaternion global radii search the complete unit similarity class',@radius; ...
 'Native delay rays and radii preserve declared grade and input coefficients',@delay; ...
 'Native varying-leading and proportional families retain regularity events',@leading; ...
 'Native quaternion state matrices use the discrete generalized Stein path',@matrix ...
};
end
function ray()
P=quaternion(reshape([-.5 1],1,1,[]),zeros(1,1,2),zeros(1,1,2),zeros(1,1,2));B=quaternion(1,0,0,0);
[~,R] = qstab.raymargin(P,B,quaternion(1,0,0,0),'Variable','z');assert(strcmp(R.status,'computed'));
assert(abs(R.negative+.5)<1e-9&&abs(R.positive-1.5)<1e-9);
end
function quotients()
P=quaternion([.1 .2 1],[.01 .04 0],[.02 -.03 0],[0 .01 0]);B=quaternion([.4 .3],[.2 .1],[.1 -.2],[-.1 .2]);
% Nominal P is close to the real Schur quadratic z^2+.2z+.1.
for side={'left','right'}
 [r,I]=qstab.freqmargin(P,B,[.2 .9 2.1],'Variable','z','Side',side{1});assert(all(I.valid)&&all(isfinite(r)));
 for j=1:numel(r),assert(I.witnesses{j}.valid&&I.witnesses{j}.evaluationResidual<1e-8);end
end
end
function radius()
P=quaternion([-.2 1],[0 0],[-.3 0],[0 0]);
[~,R] = qstab.stabilityradius(P,quaternion(1,0,0,0),'Variable','z','RelTol',.005,'AbsTol',1e-5,'MaxCells',2047);
r=1-sqrt(.13);assert(strcmp(R.status,'converged')&&R.lower<=r+1e-9&&R.upper>=r-1e-9&&R.witness.valid);
end
function delay()
P=quaternion([1 -.5],[0 0],[0 0],[0 0]);B=quaternion([0 1],[0 0],[0 0],[0 0]);savedP=P;savedB=B;
[~,R] = qstab.raymargin(P,B,1,'Variable','zinv','Grade',2);assert(strcmp(R.status,'computed')&&R.grade==2);
[~,G] = qstab.stabilityradius(P,B,'Variable','zinv','Grade',2,'RelTol',.005);assert(strcmp(G.status,'converged')&&G.lower<=.5&&G.upper>=.5);
assert(isequal(P,savedP)&&isequal(B,savedB)&&size(G.witness.coefficientComponents,2)==3);
end
function leading()
P=quaternion([-.5 1],[0 0],[0 0],[0 0]);
[~,R] = qstab.raymargin(P,P,1,'Variable','z');assert(strcmp(R.negativeEvent,'degreeLoss')&&abs(R.negative+1)<1e-12);
[~,G] = qstab.stabilityradius(P,P,'Variable','z','RelTol',.005);assert(strcmp(G.limitingEvent,'degreeLoss')&&G.lower<=1&&G.upper>=1-1e-12);
end
function matrix()
A=quaternion(.2,0,.3,0);B=quaternion(1,0,0,0);u=quaternion(0,0,0,1);
[~,R] = qstab.raymargin(A,B,u,'Type','matrix','Domain','discrete');assert(strcmp(R.status,'computed'));
assert(abs(R.positive-sqrt(.87))<1e-9&&abs(R.negative+sqrt(.87))<1e-9);
end
