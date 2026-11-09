function T=qstab_cases_papers_native()
%QSTAB_CASES_PAPERS_NATIVE The same fixed source data with MATLAB quaternion.
T={ ...
 'Native paper input preserves the qB01 left/right direction and first contacts',@sides; ...
 'Native paper catalogue retains the qB02 non-i critical box direction',@point; ...
 'Native qB03 samples and controls agree with explicit component inputs',@geometry; ...
 'Native source-linked zero-polynomial loss is separate from instability',@degree};
end
function E=go(id,native)
E=run_qstab_papers('Cases',id,'UseQuaternion',native,'Verbose',false,'MaxTime',Inf);
assert(E.ok&&E.fullyMatched);
end
function sides()
E=go('qb01_ray_sides',true);assert(E.nativeQuaternionInputs);
assert(abs(E.cases.output.left.positive-1)<1e-8);
end
function point()
E=go('qb02_edge_witness',true);assert(E.cases.output.point.witnesses{1}.valid);
end
function geometry()
A=go('qb03_parabola',true);B=go('qb03_parabola',false);
assert(isequaln(A.cases.output.values.points,B.cases.output.values.points));
assert(isequaln(A.cases.output.enclosure.regions,B.cases.output.enclosure.regions));
end
function degree()
E=go('qb03_degree_loss',true);assert(strcmp(E.cases.output.result.status,'assumptionFailed'));
end
