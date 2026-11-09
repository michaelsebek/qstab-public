function results=run_qb03_four_generator_examples()
%RUN_QB03_FOUR_GENERATOR_EXAMPLES Separate article-oriented qB03 demo.
% Complete subfamilies are checked, not four isolated vertices. No new data
% class or extra toolbox is used. Native quaternion input is accepted by the
% same interface; components make this runner executable in base MATLAB.

V=cell(1,6);
for j=1:6
    V{j}=struct('a',[1+.1*j,1],'b',[.02*j,0], ...
        'c',[.03*(-1)^j,0],'d',[.01*mod(j,3),0]);
end
args={'Reduction','four-generator','MaxTime',Inf};
results=struct();
[~,results.hurwitz] = qstab.zeroexclude(V,args{:});
[~,results.budget] = qstab.zeroexclude(V,args{:},'MaxSubfamilies',1);
S=cell(1,5);D=cell(1,5);
for j=1:5
    S{j}=struct('a',[-.2-.01*j,1],'b',[.01*j,0], ...
        'c',[.02*(-1)^j,0],'d',[.01*mod(j,3),0]);
    D{j}=struct('a',fliplr(S{j}.a),'b',fliplr(S{j}.b), ...
        'c',fliplr(S{j}.c),'d',fliplr(S{j}.d));
end
[~,results.schur] = qstab.robuststable(S,'Method','zeroexclusion',args{:},'Variable','z');
[~,results.delay] = qstab.zeroexclude(D,args{:},'Variable','zinv');
signRows=[1,1,1;1,-1,-1;-1,1,-1;-1,-1,1];
T=cell(1,5);
for j=1:4
    y=2*signRows(j,:);
    T{j}=struct('a',[-1,1,1],'b',[y(1),3*y(1),0], ...
        'c',[y(2),3*y(2),0],'d',[y(3),3*y(3),0]);
end
T{5}=[2,3,1];
[~,results.tetrahedron] = qstab.zeroexclude(T,args{:});
names=fieldnames(results);
fprintf('qB03 four-generator reduction (complete subfamilies)\n');
for j=1:numel(names)
    R=results.(names{j});
    fprintf('  %-12s : %-14s %g/%g subfamilies, %d cells\n', ...
        names{j},R.status,R.reduction.processedSubfamilies, ...
        R.reduction.totalSubfamilies,R.processedCells);
end
fprintf('The budget example is intentionally undetermined; the tetrahedron is intentionally notStable.\n');
end
