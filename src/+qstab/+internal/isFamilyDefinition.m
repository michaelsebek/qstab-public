function yes=isFamilyDefinition(F)
%ISFAMILYDEFINITION Is this a tagged public family definition?
yes=isstruct(F)&&isscalar(F)&&isfield(F,'schema')&&isequal(F.schema,'qstab.family-definition.1');
end
