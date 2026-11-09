function args=familyOptions(F,args,expected,require)
%FAMILYOPTIONS Reject conflicting metadata and pass compatible options onward.
% EXPECTED is a name/value list. REQUIRE decides which to add if omitted.
if ~qstab.internal.isFamilyDefinition(F),return;end
if mod(numel(args),2),error('qstab:family:Options','Supply complete name/value pairs.');end
for j=1:2:numel(expected)
 key=expected{j};v=expected{j+1};hit=find(cellfun(@(x)((ischar(x)&&isrow(x))||(isstring(x)&&isscalar(x)))&&strcmpi(x,key),args(1:2:end)));
 if isempty(hit)
  if any(strcmpi(key,require)),args(end+1:end+2)={key,v};end
 elseif numel(hit)>1
  error('qstab:family:OptionConflict','Duplicate %s option.',key);
 else
  ix=2*hit;w=args{ix};
  if ischar(v)
   ok=((ischar(w)&&isrow(w))||(isstring(w)&&isscalar(w)&&~ismissing(w)))&&strcmpi(w,v);
  else
   ok=isnumeric(w)&&isreal(w)&&isscalar(w)&&w==v;
  end
  if ~ok,error('qstab:family:OptionConflict','The %s option conflicts with the family definition.',key);end
  if any(strcmpi(key,require)),args{ix}=v;else,args(ix-1:ix)=[];end
 end
end
end
