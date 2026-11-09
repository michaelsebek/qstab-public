function T=qstab_cases_paperextensions_native()
T={ ...
 'Native source nearest-sphere witnesses agree with component inputs',@sphere; ...
 'Native nonlinear and ordered source fixtures retain every matched finite check',@identities; ...
 'Native rotating-factor zeros retain the two expected algebraic classes',@classes; ...
 'Native paper-figure calculations agree with their component data curves',@figures};
end
function E=go(id)
E=run_qstab_paper_extensions('Cases',id,'UseQuaternion',true,'Verbose',false);assert(E.ok&&E.fullyMatched&&E.nativeQuaternionInputs);
end
function sphere()
A=go('qb01_nearest_sphere');B=run_qstab_paper_extensions('Cases','qb01_nearest_sphere','UseQuaternion',false,'Verbose',false);
assert(isequal(A.cases.output.parameters,B.cases.output.parameters)&&max(A.cases.output.residuals)<1e-12);
end
function identities()
ids={'qb02_damping','qb02_nonlinear','qb02_correlated_marginals','qb03_ordered_product','qb03_feedback_order'};
for j=1:numel(ids),go(ids{j});end
end
function classes()
A=go('qb02_rotating_factors');r=A.cases.output.representatives;
assert(max(abs(r(:,1)-(-2+1i)))<1e-8&&max(abs(r(:,2)-(-1+1i)))<1e-8);
end
function figures()
A=qstab_paper_figures('UseQuaternion',true,'Samples',21);B=qstab_paper_figures('UseQuaternion',false,'Samples',21);
assert(A.ok&&B.ok&&A.nativeQuaternionInputs);
for j=1:numel(A.items)
 for k=1:numel(A.items(j).curves)
  a=A.items(j).curves(k);b=B.items(j).curves(k);assert(norm(a.x-b.x)+norm(a.y-b.y)+norm(a.z-b.z)<1e-9);
 end
end
end
