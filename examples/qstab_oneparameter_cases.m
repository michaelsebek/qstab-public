function C=qstab_oneparameter_cases()
%QSTAB_ONEPARAMETER_CASES Data for the scalar one-parameter article, revision 7.
%   C=qstab_oneparameter_cases returns the example catalogue without running
%   a solver. Coefficients in C.data are explicit four-component arrays;
%   they are converted to native quaternion inputs by run_qstab_oneparameter.
%   C.rays contains the rational directions, unchanged in length or side.
%   Recorded analytical targets are references, not newly computed answers.
%
%   E1--E7 are the seven archived frequency datasets. E5V is the additional
%   noncommuting tangency, E5C the real cubic multiplicity diagnostic, and
%   D11--D13 the three remaining directional diagnostics. A01 checks zero
%   perturbations, nonmonic data and nonunit directions.
%   See also run_qstab_oneparameter, plot_qstab_oneparameter.
root=fileparts(fileparts(mfilename('fullpath')));
F=jsondecode(fileread(fullfile(root,'reproducibility','qb01','inputs.json')));
ids={'E1','E2','E3','E4','E5','E6','E7','E5V','E5C','D11','D12','D13','A01'};
titles={'Sensor and actuator calibration', ...
 'Rank-one perturbation frequency matrix','Sphere of nearest gains', ...
 'Nearly repeated singular values','Tangential first contact', ...
 'Nearly flat frequency optimum','Amplitude and frequency rescaling', ...
 'Tangency without a common complex slice','Real cubic and inflated multiplicity', ...
 'Infinite directional margin','Later guardian candidate without imaginary contact', ...
 'Spherical contact with zero odd part','Zero, nonmonic and nonunit-direction checks'};
C=struct('id',{},'title',{},'sourceRevision',{},'sourceAnchor',{},'data',{},'rays',{});
for k=1:numel(ids)
 d=struct();ii=find(strcmp({F.datasets.id},ids{k}));
 if ~isempty(ii),d=F.datasets(ii);end
 rr=F.rays(strcmp({F.rays.articleCase},ids{k}));
 anchor=['Example ' ids{k}];
 if ids{k}(1)=='D',anchor=['reproducibility/directional/' rr.name];end
 if strcmp(ids{k},'A01'),anchor='Assumption (A); ray normalization and nonmonic invariance';end
 C(end+1)=struct('id',ids{k},'title',titles{k},'sourceRevision',F.sourceRevision, ...
  'sourceAnchor',anchor,'data',d,'rays',rr); %#ok<AGROW>
end
end
