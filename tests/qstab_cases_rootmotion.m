function cases=qstab_cases_rootmotion()
%QSTAB_CASES_ROOTMOTION Analytical simple-root jets and failure guards.
cases={'OOO contact gives the article Hessian without finite differences',@ooo; ...
 'EOO contact gives the article Hessian without finite differences',@eoo; ...
 'Repeated and nonroot representatives remain unresolved',@unresolved; ...
 'Real parameter units and common coefficient scale are respected',@units; ...
 'Out-of-box and simplex coordinates are rejected',@invalid};
end
function [P,G]=data(which)
if which==1
 C=[-1/8,-1/10,1/4,1/10,1;0,5/8,1/5,5/4,0;3/40,0,3/40,0,0;zeros(1,5)];
 powers=[1 1 1];
else
 C=[-1,-1/2,-1/4,3/2,1;0,15/8,2,17/8,0;0,-3/2,0,-3/2,0;zeros(1,5)];
 powers=[0 1 1];
end
P=parts(C);G=cell(1,3);
for k=1:3,D=zeros(4,5);D(k+1,powers(k)+1)=1;G{k}=parts(D);end
end
function ooo()
[P,G]=data(1);F=qstab.family('affine',P,G,repmat([-.001,.001],3,1));
[g,H,I]=qstab.rootmotion(F,1i);
Q=(64/93015)*[2612,-1125,0;-1125,900,0;0,0,900];
assert(I.valid&&norm(g-[-8i/3;0;0])<1e-11&&norm(-real(H)/2-Q,'fro')<1e-10);
end
function eoo()
[P,G]=data(2);F=qstab.family('affine',P,G,repmat([-.001,.001],3,1));
[g,H,I]=qstab.rootmotion(F,1i);
Q=[209,-48,0;-48,16,0;0,0,16]/520;
assert(I.valid&&norm(g-[-1i;0;0])<1e-11&&norm(-real(H)/2-Q,'fro')<1e-10);
end
function unresolved()
F=qstab.family('affine',[1,2,1],{1},[-.1,.1]);
[g,H,I]=qstab.rootmotion(F,-1);assert(~I.valid&&all(isnan([g(:);H(:)])));
[~,~,I]=qstab.rootmotion(F,2);assert(~I.valid);
end
function units()
[P,G]=data(1);F=qstab.family('affine',P,G,repmat([-.001,.001],3,1));
[g,H]=qstab.rootmotion(F,1i);C=qstab.polycomponents(P)*1e-100;
S=diag([2,3,4]);GG=G;
for k=1:3,GG{k}=parts(qstab.polycomponents(G{k})*S(k,k)*1e-100);end
F=qstab.family('affine',parts(C),GG,repmat([-.0001,.0001],3,1));
[gs,Hs,I]=qstab.rootmotion(F,1i);
assert(I.valid&&norm(gs-S*g)<1e-10&&norm(Hs-S*H*S,'fro')<1e-8);
end
function invalid()
[P,G]=data(1);F=qstab.family('affine',P,G,repmat([-.001,.001],3,1));
caught=false;try,qstab.rootmotion(F,1i,'Parameters',[1;0;0]);catch e,caught=strcmp(e.identifier,'qstab:rootmotion:Parameters');end;assert(caught);
caught=false;try,qstab.rootmotion(qstab.family('polytope',{P,P}),1i);catch e,caught=strcmp(e.identifier,'qstab:rootmotion:Family');end;assert(caught);
end
function P=parts(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
