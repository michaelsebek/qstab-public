function cases=qstab_cases_intervalstudy()
%QSTAB_CASES_INTERVALSTUDY Reproduce current article fixtures and run contracts.
cases={'All 24 current interval article/supplement groups compute their references',@allCases; ...
 'Selected interval examples preserve path, folder and RNG',@selected; ...
 'Invalid selections and unknown options fail before any example',@invalid};
end
function allCases()
R=run_qstab_intervals('UseQuaternion',false,'Save',false,'Verbose',false);
assert(R.ok&&numel(R.cases)==24&&R.checksPassed>=95&&~R.exact.executed);
assert(all(strcmp({R.cases.status},'MATCHED')));
end
function selected()
p=path;d=pwd;s=rng;
R=run_qstab_intervals({'I09','I20'},'UseQuaternion',false,'Save',false,'Verbose',false);
assert(R.ok&&numel(R.cases)==2&&strcmp(path,p)&&strcmp(pwd,d)&&isequaln(rng,s));
assert(~R.exact.executed&&~R.exact.ok);
end
function invalid()
p=path;d=pwd;
caught=false;try,run_qstab_intervals({'I01','I01'},'Save',false);catch e,caught=strcmp(e.identifier,'qstab:intervalStudy:Selection');end;assert(caught);
caught=false;try,run_qstab_intervals('I01','Unknown',1);catch,caught=true;end;assert(caught&&strcmp(path,p)&&strcmp(pwd,d));
end
