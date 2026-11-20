function [D,checks]=qstab_interval_evaluate(id,o)
%QSTAB_INTERVAL_EVALUATE Calculations, independent references, and evidence.
% Called by run_qstab_intervals. References are never used as solver outputs.
checks=struct('name',{},'actual',{},'reference',{},'tolerance',{},'pass',{});
D=struct();
switch id
case 'I01'
 C=[1 1;1 0;0 0;0 0];P=pol(C,o);
 D.norm=qstab.normpoly(P);[D.representatives,D.spectrum]=qstab.roots(P);
 D.actualRoot=[-1;-1;0;0];D.evaluation=evalq(C,D.actualRoot);
 [D.stable,D.stability]=qstab.isstable(P);
 add('Real central norm',D.norm,[2 2 1],1e-13);add('Canonical class',D.representatives,-1+1i,1e-10);
 add('Actual right zero residual',norm(D.evaluation),0,1e-13);add('Nominal stability',D.stable,1,0);
case 'I02'
 D.parameters=linspace(1,2,21);D.normCoefficients=zeros(21,3);
 for k=1:21,C=[D.parameters(k),1;1,0;zeros(2,2)];D.normCoefficients(k,:)=qstab.normpoly(pol(C,o));end
 D.couplingResidual=D.normCoefficients(:,1)-1-(D.normCoefficients(:,2)/2).^2;
 add('Quadratic coupling is preserved',max(abs(D.couplingResidual)),0,1e-13);
case 'I03'
 C=[1 2 1;zeros(3,3)];W=zeros(4,3);W(2,2)=1;W(3,1)=1;
 [D.radius,D.radiusGauge]=qstab.boxvalue(pol(C,o),W,points([0;0;-1;0],o),'Witness',true);
 D.rhos=[1.5 2 3];D.normValueRanges=zeros(3,2);D.cornerFlags=zeros(4,1);
 for k=1:3,r=D.rhos(k);D.normValueRanges(k,:)=[-4-r^2,-4+r^2];end
 D.frequency=linspace(0,3,301);D.frequencyRanges=cell(3,1);
 for k=1:3,r=D.rhos(k);w=D.frequency;D.frequencyRanges{k}=[(w.^4-(6+r^2)*w.^2+1).', (w.^4-6*w.^2+1+r^2).', (4*w.*(1-w.^2)).'];end
 params=[-3 -3;-3 3;3 -3;3 3];D.cornerNorms=zeros(4,5);
 for k=1:4,Q=C;Q(2,2)=params(k,1);Q(3,1)=params(k,2);D.cornerFlags(k)=qstab.isstable(pol(Q,o));D.cornerNorms(k,:)=qstab.normpoly(pol(Q,o));end
 Q=C;Q(3,1)=3;D.unstable=nominal(Q,'s',o);D.unstableRightRoot=[-1+sqrt(3/2);0;-sqrt(3/2);0];
 add('All four vertices stable',D.cornerFlags,ones(4,1),0);add('Edge member unstable',D.unstable.flag,0,0);
 add('Edge right zero',norm(evalq(Q,D.unstableRightRoot)),0,1e-12);
 Q=C;Q(3,1)=2;D.contact=nominal(Q,'s',o);D.contactRightRoot=[0;0;-1;0];
 D.nonzeroAtI=evalq(Q,[0;1;0;0]);add('Contact factorization norm',D.contact.norm,[5 4 6 4 1],1e-12);
 add('Contact right zero',norm(evalq(Q,D.contactRightRoot)),0,1e-12);add('Representative need not evaluate to zero',D.nonzeroAtI,[0;2;2;0],1e-12);
 D.edgeContacts=[-sqrt(5),sqrt(5)];D.edgeContactResidual=zeros(1,2);
 for k=1:2,Q=C;Q(2,2)=D.edgeContacts(k);Q(3,1)=3;u=[0;0;-2/3;D.edgeContacts(k)/3];D.edgeContactResidual(k)=norm(evalq(Q,u));end
 add('Both continuous-edge contacts',D.edgeContactResidual,[0 0],1e-12);
 F=qstab.family('interval',pol(C-3*W,o),pol(C+3*W,o));[D.cube,D.cubeInfo]=qstab.kharitonovcube(F);
 add('Two active cube coordinates',D.cubeInfo.activeParameterCount,2,0);
 [D.independentFlag,D.independentReport]=qstab.robuststable(F,'Method','lyapunov','MaxCells',256,'MaxTime',10);
 add('Running-square first contact radius',D.radius,2,1e-10);
 add('Independent square is not robustly stable',D.independentFlag,0,0);
 G=zeros(4,3);G(2,2)=1;G(3,1)=1;D.correlated=qstab.family('affine',pol(C,o),{pol(G,o)},[-3,3]);
 D.correlatedParameters=[-3 -1 0 1 3];D.correlatedFlags=zeros(size(D.correlatedParameters));
 for k=1:numel(D.correlatedParameters),D.correlatedFlags(k)=qstab.isstable(pol(C+D.correlatedParameters(k)*G,o));end
 add('Correlated members',D.correlatedFlags,ones(size(D.correlatedFlags)),0);
 add('Source norm-value ranges',D.normValueRanges,[-6.25,-1.75;-8,0;-13,5],1e-12);
 D.continuumEvidence='The formula beta^2<4+alpha^2 proves the claims; samples alone do not.';
case 'I04'
 D.parameters=[0 1 2];D.members=cell(1,3);
 for k=1:3,D.members{k}=nominal([1,1-D.parameters(k);zeros(3,2)],'s',o,1);end
 add('Stable before degree loss',D.members{1}.flag,1,0);add('Unstable after degree loss',D.members{3}.flag,0,0);
 add('Infinity retained at declared grade',sum(isinf(D.members{2}.roots)),1,0);
case 'I05'
 D.alpha=[9/20,3/4,3/2];D.flags=zeros(1,3);D.delta3=zeros(1,3);D.members=cell(1,3);
 for k=1:3,a=D.alpha(k);C=[2,a,1;4,0,0;1,1,0;0,2,0];D.members{k}=nominal(C,'s',o);D.flags(k)=D.members{k}.flag;D.delta3(k)=D.members{k}.hurwitzDeterminants(3);end
 add('Stable endpoints and unstable interior',D.flags,[1 0 1],0);add('Third Hurwitz determinants',D.delta3,[451/20000,-17/32,17],5e-10);
case 'I06'
 C=[1 1 1;2 0 0;zeros(2,3)];D.centre=nominal(C,'s',o);
 D.boundaryY=[4 6 8];D.boundaryMembers=cell(1,3);D.boundaryFlags=zeros(1,3);
 for k=1:3,Q=C;Q(3,2)=2;Q(4,2)=sqrt(D.boundaryY(k)-4);D.boundaryMembers{k}=nominal(Q,'s',o);D.boundaryFlags(k)=D.boundaryMembers{k}.flag;end
 D.boundaryMinimumDelta3=4*(4-3);D.wholeValueRange=[-5 3];D.boundaryValueRange=[-5 -1];
 add('Unstable centre',D.centre.flag,0,0);add('Centre Routh first column',D.centre.routhColumn,[1 2 2 -3 5].',1e-12);
 add('Boundary members',D.boundaryFlags,[1 1 1],0);add('Whole-boundary determinant lower bound',D.boundaryMinimumDelta3,4,0);
case 'I07'
 C=[1 2 1;0 3 0;3 0 0;zeros(1,3)];W=zeros(4,3);W(2,2)=1;
 [D.radius,D.frequencyReport]=qstab.boxfreqmargin(pol(C,o),W,1,'Witness',true);
 D.reference=3-sqrt(5);add('Slice frequency/contact reserve',D.radius,D.reference,1e-10);
 D.frequency=linspace(.05,2,300);D.frequencyMargin=qstab.boxfreqmargin(pol(C,o),W,D.frequency);
 D.globalReason='The only positive contact frequency is omega=1; this independent slice has beta fixed at 3.';
case 'I08'
 C=[1 2 1;zeros(3,3)];W=zeros(4,3);W(2,2)=1;W(3,1)=1;
 D.directions=[0,0,0,0,0,0;0,0,0,1,0,3/5;1,-1,3/5,0,0,4/5;0,0,4/5,0,1,0];
 [D.radius,D.gauge]=qstab.boxvalue(pol(C,o),W,points(D.directions,o),'Witness',true);
 add('Point gauges, including infeasibility',D.radius,[2,2,10/3,Inf,Inf,Inf],1e-9);
case {'I09','I10'}
 which=1+strcmp(id,'I10');[C,G,delta,Qref,gref]=quartic(which);
 generators=cell(1,3);for k=1:3,generators{k}=pol(G(:,:,k),o);end
 F=qstab.family('affine',pol(C,o),generators,repmat([-.001,.001],3,1));
 [D.gradient,D.Hessian,D.motion]=qstab.rootmotion(F,1i);
 D.quadraticForm=-real(D.Hessian)/2;D.base=nominal(C,'s',o);D.shift=delta;
 D.shiftedCoefficients=shiftpoly(C,-delta);D.shiftedCentre=nominal(D.shiftedCoefficients,'s',o);
 D.boundaryParameters=boundarySamples(3,.001,3);D.boundaryAbscissae=zeros(size(D.boundaryParameters,1),1);
 for k=1:size(D.boundaryParameters,1),T=D.shiftedCoefficients;for j=1:3,T=T+D.boundaryParameters(k,j)*G(:,:,j);end;rr=qstab.roots(pol(T,o));D.boundaryAbscissae(k)=max(real(rr));end
 [D.faceMinima,D.minimizers]=quadraticBoxMin(D.quadraticForm);
 D.boundaryQuadraticMinimum=min(D.faceMinima);D.exactProofKey={'quarticOOO','quarticEOO'};D.exactProofKey=D.exactProofKey{which};
 D.baseFactorNorm=D.base.norm;D.simpleContactNormDerivative=D.motion.normDerivative*D.motion.coefficientScale^2;
 add('Simple root is verified',D.motion.valid,true,0);add('Analytical first derivatives',D.gradient,gref,1e-10);
 add('Contact quadratic form',D.quadraticForm,Qref,1e-10);
 refs=[560/1959,2/209];add('Minimum of quadratic form on cube boundary',D.boundaryQuadraticMinimum,refs(which),1e-11);
 add('Shifted central rightmost root',max(real(D.shiftedCentre.roots)),delta,2e-11);
 D.boundaryEvidence='Displayed roots are samples. Full boundary proof is replayed separately with exact rational Bernstein arithmetic.';
case {'I11','I12','I13','I14'}
 kind=str2double(id(2:end))-10;which=1+mod(kind-1,2);[C,G,delta]=quartic(which);C=shiftpoly(C,-delta);epsilon=1e-12;
 if kind<=2,C(:,6)=0;C(1,6)=epsilon;G(:,6,:)=0;
 else,C=[zeros(4,1),C];C(1,1)=-epsilon;G=cat(2,zeros(4,1,3),G);end
 D.coefficients=C;D.epsilon=epsilon;D.centre=nominal(C,'s',o,5);D.normDegree=numel(D.centre.norm)-1;
 D.normalizedRouthColumn=routh(D.centre.norm/max(abs(D.centre.norm)));
 D.negativeIndex=8+2*(kind<=2);D.normalizedNegativeEntry=D.normalizedRouthColumn(D.negativeIndex);
 col=D.centre.routhColumn;D.floatingRouthSignChanges=sum(sign(col(1:end-1)).*sign(col(2:end))<0);
 D.negativeEntry=col(D.negativeIndex);
 D.floatingRouthResolved=all(isfinite(col))&&all(col~=0);
 D.boundaryParameters=boundarySamples(3,.001,2);D.sampledRoots=cell(size(D.boundaryParameters,1),1);
 for k=1:size(D.boundaryParameters,1),T=C;for j=1:3,T=T+D.boundaryParameters(k,j)*G(:,:,j);end;D.sampledRoots{k}=qstab.roots(pol(T,o));end
 names={'OOO','EOO','EEE','EEO'};D.kind=names{kind};D.exactProofKey=['quintic' D.kind];
 add('No tiny coefficient is trimmed',D.normDegree,10,0);add('Declared grade remains five',D.centre.rootInfo.grade,5,0);
 add('Every floating Routh entry is finite and nonzero',D.floatingRouthResolved,true,0);
 add('Floating Routh sign changes agree with the exact reference',D.floatingRouthSignChanges,2,0);
 if which==1,bounds=[-1/500000,-1/1000000];else,bounds=[-7/100000000,-3/50000000];end
 D.paperNegativeEntryBounds=bounds;
 add('Floating negative entry is inside the paper rational bounds',D.negativeEntry>bounds(1)&&D.negativeEntry<bounds(2),true,0);
 D.proofNotice='Exact centre Routh counts and all six continuous boundary faces are verified by verify_qstab_interval_certificates / reproducibility/qb02. Default floating-point flags remain unchanged.';
case 'I15'
 C=[1 4 6 4 1;zeros(3,5)];W=zeros(4,5);W(1,1)=1/10;W(2,2)=1;
 U=[0,0,0,0;1,-1,0,3/5;0,0,1,4/5;0,0,0,0];
 [D.radii,D.gauges]=qstab.boxvalue(pol(C,o),W,points(U,o),'Witness',true);
 [D.zeroRadius,D.zeroReport]=qstab.boxvalue(pol(C,o),W,0);
 D.ranks=zeros(1,4);for k=1:4,D.ranks(k)=D.gauges.pointData{k}.gauge.numericalRank;end
 Q=C;Q(1,1)=Q(1,1)+4/11;Q(2,2)=-40/11;D.contact=nominal(Q,'s',o);D.rightRoot=[0;1;0;0];
 add('Exceptional versus regular directions',D.radii,[40/11,40/11,40,40],1e-8);
 add('Generator rank',D.ranks,[1 1 2 2],0);add('Zero-frequency radius',D.zeroRadius,10,1e-9);
 add('Actual contact right zero',norm(evalq(Q,D.rightRoot)),0,1e-12);
 D.globalRadius=min([D.zeroRadius,D.radii]);D.globalReason='Coefficient-slice imaginary part 4*w*(1-w^2) leaves only w=0 and w=1; compare their gauges.';
case 'I16'
 D.degrees=[1 2 5];D.members=cell(3,2);
 for k=1:3,n=D.degrees(k);C=zeros(4,n+1);C(1,n)=1.5;C(1,n+1)=1;D.members{k,1}=nominal(C,'z',o,n);C(2,n+1)=2;D.members{k,2}=nominal(C,'z',o,n);
  add(sprintf('Degree %d centre unstable',n),D.members{k,1}.flag,0,0);add(sprintf('Degree %d boundary stable',n),D.members{k,2}.flag,1,0);
  add(sprintf('Degree %d zero modes',n),sum(abs(D.members{k,2}.roots)<1e-10),n-1,0);
 end
 D.boundaryRootModulusBound=1.5/sqrt(5);add('Entire cube-boundary modulus bound below one',D.boundaryRootModulusBound<1,true,0);
case 'I17'
 C=[0 0 1;0 6/5 0;zeros(2,3)];D.centre=nominal(C,'z',o);
 D.boundarySquaredNorm=[49/100,98/100];D.boundary=cell(1,2);
 for k=1:2,Q=C;Q(3,1)=.7;Q(4,1)=sqrt(D.boundarySquaredNorm(k)-.49);D.boundary{k}=nominal(Q,'z',o);add('Boundary Schur member',D.boundary{k}.flag,1,0);end
 D.completeBoundaryGaps=[min(D.boundarySquaredNorm)-11/25,1-max(D.boundarySquaredNorm)];
 add('Centre is unstable',D.centre.flag,0,0);add('Whole-boundary annulus gaps',D.completeBoundaryGaps,[1/20,1/50],1e-14);
 Q=C;Q(3,1)=sqrt(11/25);D.innerContactNorm=qstab.normpoly(pol(Q,o));add('Inner contact factorization',D.innerContactNorm,[11/25,0,36/25,0,1],1e-12);
case 'I18'
 C=[-13/50,-3/4,0,1;zeros(3,4)];D.centre=nominal(C,'z',o);D.v=[1/4,3/4];D.boundary=cell(1,2);
 B=-3/4;cc=-13/50;dd=1-cc^2;uu=-B*cc;vv=B;
 D.L=zeros(1,2);D.F=zeros(1,2);
 for k=1:2,x=D.v(k);Q=C;Q(2,3)=sqrt(x);D.boundary{k}=nominal(Q,'z',o);D.L(k)=dd^2-vv^2-cc^2*x;D.F(k)=D.L(k)^2-(uu*(dd-vv)-cc*x)^2-x*(dd+vv-cc*uu)^2;add('Schur cubic boundary member',D.boundary{k}.flag,1,0);end
 D.FquadraticCoefficient=cc^4-cc^2;D.centreAtOne=sum(C(1,:));
 add('Unstable cubic centre',D.centre.flag,0,0);add('Centre p(1)',D.centreAtOne,-1/100,1e-15);
 add('Complete radial interval L minimum',D.L(2),1601061/6250000,1e-14);
 add('Concave F endpoint values',D.F,[1647104154759,5454247783509]/156250000000000,1e-13);
 add('Concavity of F',D.FquadraticCoefficient<0,true,0);
case 'I19'
 D.marginPolynomial=[51 -920 5500 -5000];D.radius=fzero(@(r)polyval(D.marginPolynomial,r),[1 1.1]);
 D.isolatingInterval=[1098721768495,1098721768496]/1e12;D.gapAtBudgets=polyval([-51/1000,23/25,-11/2,5],[1 1.1]);
 r=D.radius;d=1-r/10;kap=d;b=2+r/10;y=2-r/2;z=y;S=d*d+y*y+z*z;
 Q=[kap,d,1;b,0,0;0,y,0;0,z,0];D.contact=nominal(Q,'s',o);D.rightRoot=(b/S)*[0;-d;z;-y];D.physicalContactFrequency=5*sqrt(kap);
 D.budgets=[0 1 1.2];D.time=(0:.1:40).';D.trajectories=zeros(numel(D.time),3);D.spectralAbscissae=zeros(1,3);D.matrices=cell(1,3);
 for k=1:3,r=D.budgets(k);d=1-r/10;kap=d;b=2+r/10;y=2-r/2;z=y;
  A=[zeros(4),eye(4);-25*leftmatrix([kap;b;0;0]),-5*leftmatrix([d;0;y;z])];D.matrices{k}=A;D.spectralAbscissae(k)=max(real(eig(A)));
  for j=1:numel(D.time),x=expm(A*D.time(j))*[10;zeros(7,1)];D.trajectories(j,k)=norm(x(1:4));end
 end
 add('Five-gain reserve',D.radius,1.0987217684952284,2e-13);add('Gap at budgets 1 and 1.1',D.gapAtBudgets,[369/1000,-4681/1e6],2e-14);
 add('Actual quaternion root',norm(evalq(Q,D.rightRoot)),0,1e-12);
 add('Physical spectral abscissae',D.spectralAbscissae,[-.2620761414007129,-.0524051008658224,.0658486838669244],1e-11);
 D.physicalGainsAtBudgetOne=[.9,1.1;4.5,5.5;9.5,10.5;1.5,2.5;1.5,2.5];
case 'I20'
 C=[1 1 1;2 0 0;0 2 0;0 2 0];W=zeros(4,3);W(1,1:2)=.1;
 D.u=fzero(@(u)u^3+8*u-4,[4/9,.5]);D.radius=10*(1-D.u);D.omega=sqrt(D.u);
 D.frequency=sort([linspace(.15,1.2,400),D.omega,sqrt(.5)]);
 [D.frequencyMargin,D.frequencyInfo]=qstab.boxfreqmargin(pol(C,o),W,D.frequency);
 [D.contactRadius,D.contactReport]=qstab.boxfreqmargin(pol(C,o),W,D.omega,'Witness',true);
 D.zeroFrequencyRadius=qstab.boxfreqmargin(pol(C,o),W,0);
 u=D.u;Q=C;Q(1,1:2)=u;D.contact=nominal(Q,'s',o);D.rightRoot=[0;-u*u/2;u;-u];
 D.gapAtOne=.9*(.9^2+8)-4;D.reservePolynomialResidual=polyval([1 -30 1100 -5000],D.radius);
 add('Restricted reserve',D.radius,5.143205348290382,1e-12);add('Frequency optimum',D.contactRadius,D.radius,1e-10);
 add('Critical frequency',D.omega,.6969070706851537,1e-13);add('Zero frequency unattainable',D.zeroFrequencyRadius,Inf,0);
 add('Right root residual',norm(evalq(Q,D.rightRoot)),0,1e-12);add('Budget-one gap',D.gapAtOne,3929/1000,1e-13);
 add('Algebraic reserve polynomial',D.reservePolynomialResidual,0,2e-10);
case 'I21'
 C=[1 1 1;2 0 0;0 2 0;zeros(1,3)];W=zeros(4,3);W(3,2)=1;W(4,2)=2;
 D.rightRoot=[0;-.5;0;-sqrt(3)/2];[D.radius,D.gauge]=qstab.boxvalue(pol(C,o),W,points(D.rightRoot,o),'Witness',true);
 D.vertexRayParameters=[-1 -1;-1 1;1 -1;1 1];
 D.vertexRayMinimizingRadii=max(0,-4*D.vertexRayParameters(:,1)/10);
 D.vertexRayMinima=(2+D.vertexRayParameters(:,1).*D.vertexRayMinimizingRadii).^2+4*D.vertexRayMinimizingRadii.^2;
 D.vertexRayMinimum=min(D.vertexRayMinima);D.reservePolynomialResidual=polyval([1 -4 1],D.radius);
 add('Edge-interior reserve',D.radius,2-sqrt(3),1e-10);add('All vertex rays stay outside unstable disk',D.vertexRayMinimum>3,true,0);
 add('Algebraic reserve relation',D.reservePolynomialResidual,0,1e-10);
case 'I22'
 C=[1 1;zeros(3,2)];W=zeros(4,2);W(1,1)=1;
 [D.zeroRadius,D.zeroWitness]=qstab.boxvalue(pol(C,o),W,0,'Witness',true);
 D.leadingLoss=1;D.leadingLossMember=nominal([1 0;zeros(3,2)],'s',o,1);
 D.imaginaryParameters=[-100,-2,0,2,100];D.imaginaryRealParts=zeros(1,5);
 for k=1:5,Q=C;Q(2,1)=D.imaginaryParameters(k);r=qstab.roots(pol(Q,o));D.imaginaryRealParts(k)=real(r);end
 add('Zero-root reserve',D.zeroRadius,1,1e-12);add('Degree-loss infinity',sum(isinf(D.leadingLossMember.roots)),1,0);
 add('Imaginary constant perturbations preserve real parts',D.imaginaryRealParts,-ones(1,5),1e-12);
 D.unboundedReserve=Inf;D.unboundedReason='All actual roots are -1-t*i for real t, so their real part is identically -1.';
case 'I23'
 C=[1 1;zeros(3,2)];W=.1*ones(4,2);F=qstab.family('interval',pol(C-W,o),pol(C+W,o));
 [D.linearFlag,D.linearReport]=qstab.robuststable(F,'Method','zeroexclude','MaxCells',2048,'MaxTime',15);
 D.innerProductLowerBound=.9^2-3*.1^2;add('Analytic linear box lower bound',D.innerProductLowerBound,39/50,1e-14);
 C=[1 2 1;zeros(3,3)];W=[.01*ones(4,2),zeros(4,1)];F=qstab.family('interval',pol(C-W,o),pol(C+W,o));
 [D.quadraticFlag,D.quadraticReport]=qstab.robuststable(F,'Method','zeroexclude','MaxCells',2048,'MaxTime',15);
 D.discriminant=(-2)^2-4*3;add('Quadratic strict frequency domination',D.discriminant,-8,0);
 add('Numerical linear family test',D.linearFlag,1,0);add('Numerical quadratic family test',D.quadraticFlag,1,0);
case 'I24'
 D.degrees=1:5;D.counts=zeros(1,5);D.dimensions=zeros(1,5);
 for n=1:5
  C=zeros(4,n+1);C(1,:)=fliplr(poly(-ones(1,n)));W=.001*ones(4,n+1);W(:,end)=0;
  if n==1,W(:,end)=.001;end
  F=qstab.family('interval',pol(C-W,o),pol(C+W,o));[faces,I]=qstab.intervalfaces(F);D.counts(n)=numel(faces);D.dimensions(n)=I.parameterDimensionBound;
 end
 D.fullCubeCounts=[2^8,8*2^7,nchoosek(4,2)*4*2^6,nchoosek(4,3)*8*2^5];
 D.quarticAlgebraicCandidateCount=256+5120*4+512*(4-1);
 add('Degree-dependent continuous family counts',D.counts,[256 1536 1536 896 1024],0);
 add('Full cube counts: vertices edges two-faces three-faces',D.fullCubeCounts,[256 1024 1536 1024],0);
 add('Quartic refined count',D.counts(4),896,0);add('Quartic total algebraic check count',D.quarticAlgebraicCandidateCount,22272,0);
 % The degree-independent hierarchy above and quartic-specific pruning are
 % different counts; the article's optimized quartic number is 21248.
 D.quarticNegativeThreeFaces=(nchoosek(4,3)+4*nchoosek(3,2))*2^5;
 D.optimizedQuarticCandidateCount=D.fullCubeCounts(1)+ ...
     2*4*(D.fullCubeCounts(2)+D.fullCubeCounts(3))+D.quarticNegativeThreeFaces;
 add('Optimized quartic nominal/candidate count',D.optimizedQuarticCandidateCount,21248,0);
otherwise
 error('qstab:intervalStudy:Case','Unknown interval case %s.',id);
end
    function add(name,actual,reference,tol)
        a=actual;r=reference;
        pass=isequal(size(a),size(r));
        if pass
            finite=isfinite(a)&isfinite(r);special=~finite;
            pass=all(abs(a(finite)-r(finite))<=tol.*max(1,abs(r(finite)))) && ...
                 all((isinf(a(special))&isinf(r(special))&sign(a(special))==sign(r(special))) | (isnan(a(special))&isnan(r(special))));
        end
        checks(end+1)=struct('name',name,'actual',a,'reference',r,'tolerance',tol,'pass',logical(pass));
    end
end
function P=pol(C,o)
if o.UseQuaternion,P=quaternion(C(1,:),C(2,:),C(3,:),C(4,:));
else,P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));end
end
function P=points(C,o)
P=pol(C,o);
end
function R=nominal(C,v,o,grade)
if nargin<4,grade=size(C,2)-1;end
P=pol(C,o);[R.flag,R.stability]=qstab.isstable(P,'Variable',v,'Grade',grade);
[R.roots,R.rootInfo]=qstab.roots(P,'Variable',v,'Grade',grade);
R.norm=qstab.normpoly(P,'Grade',grade);R.coefficients=C;
R.hurwitzDeterminants=hurwitz(R.norm);R.routhColumn=routh(R.norm);
end
function D=hurwitz(c)
n=numel(c)-1;a=fliplr(c);H=zeros(n);D=zeros(1,n);
for j=1:n,for k=1:n,r=2*k-j;if r>=0&&r<=n,H(j,k)=a(r+1);end;end;end
for j=1:n,D(j)=det(H(1:j,1:j));end
end
function col=routh(c)
a=fliplr(c);n=numel(a)-1;m=ceil((n+1)/2);T=zeros(n+1,m);
T(1,1:numel(a(1:2:end)))=a(1:2:end);
if n==0,col=T(1);return;end
T(2,1:numel(a(2:2:end)))=a(2:2:end);
for j=3:n+1
 if T(j-1,1)==0,T(j:end,:)=NaN;break;end
 for k=1:m-1,T(j,k)=(T(j-1,1)*T(j-2,k+1)-T(j-2,1)*T(j-1,k+1))/T(j-1,1);end
end
col=T(:,1);
end
function y=evalq(C,q)
y=C(:,end);for k=size(C,2)-1:-1:1,y=mul(y,q)+C(:,k);end
end
function c=mul(a,b)
c=[a(1)*b(1)-a(2:4).'*b(2:4);a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
end
function A=leftmatrix(q)
a=q(1);b=q(2);c=q(3);d=q(4);A=[a,-b,-c,-d;b,a,-d,c;c,d,a,-b;d,-c,b,a];
end
function [C,G,delta,Q,g]=quartic(which)
if which==1
 C=[-1/8,-1/10,1/4,1/10,1;0,5/8,1/5,5/4,0;3/40,0,3/40,0,0;zeros(1,5)];p=[1,1,1];delta=1/4000000;
 Q=64/93015*[2612,-1125,0;-1125,900,0;0,0,900];g=[-8i/3;0;0];
else
 C=[-1,-1/2,-1/4,3/2,1;0,15/8,2,17/8,0;0,-3/2,0,-3/2,0;zeros(1,5)];p=[0,1,1];delta=1/125000000;
 Q=[209,-48,0;-48,16,0;0,0,16]/520;g=[-1i;0;0];
end
G=zeros(4,5,3);for j=1:3,G(j+1,p(j)+1,j)=1;end
end
function D=shiftpoly(C,a)
n=size(C,2)-1;D=zeros(size(C));
for k=0:n,for j=0:k,D(:,j+1)=D(:,j+1)+C(:,k+1)*nchoosek(k,j)*a^(k-j);end;end
end
function B=boundarySamples(m,r,n)
% Finite plot/inspection points only, never used as a continuum proof.
v=linspace(-r,r,n);A=cell(1,m);[A{:}]=ndgrid(v);B=zeros(numel(A{1}),m);
for j=1:m,B(:,j)=A{j}(:);end
B=B(any(abs(B)==r,2),:);
end
function [values,minimizers]=quadraticBoxMin(Q)
% Exhaust all active sets on each positive face of a positive definite
% three-dimensional quadratic. Central symmetry supplies the negative faces.
m=size(Q,1);values=inf(1,m);minimizers=zeros(m);
for j=1:m
 k=setdiff(1:m,j);
 for code=0:3^(m-1)-1
  state=mod(floor(code./3.^(0:m-2)),3);free=k(state==0);fixed=k(state~=0);
  z=zeros(m,1);z(j)=1;z(fixed)=2*state(state~=0)-3;
  given=[j,fixed];if ~isempty(free),z(free)=-Q(free,free)\(Q(free,given)*z(given));end
  if any(abs(z)>1+1e-12),continue;end
  value=z.'*Q*z;if value<values(j),values(j)=value;minimizers(:,j)=z;end
 end
end
end
