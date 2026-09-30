function require_columns(T, columns, source)
%REQUIRE_COLUMNS Error if table T lacks any of the named columns.
names = T.Properties.VariableNames;
missing = columns(~ismember(columns, names));
if ~isempty(missing)
    error('focal_stress:missingColumns', '%s: missing column(s) %s. Found: %s', ...
        source, strjoin(missing, ', '), strjoin(names, ', '));
end
end
