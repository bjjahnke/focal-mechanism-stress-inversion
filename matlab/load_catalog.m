function T = load_catalog(path_or_table)
%LOAD_CATALOG Read and validate a focal mechanism catalog.
%   T = load_catalog('catalog.csv')  reads a CSV file.
%   T = load_catalog(T)              validates a table already in memory.
%
%   Required columns : strike, dip, rake (degrees)
%   Optional columns : event_id (created as EV001, EV002, ... if missing)
%                      uncertainty_deg (enables the variance test in step 3)
%   Any other columns are kept unchanged.
if istable(path_or_table)
    T = path_or_table; source = 'catalog';
else
    source = char(path_or_table);
    if ~isfile(source)
        error('focal_stress:fileNotFound', 'Catalog file not found: %s', source);
    end
    T = readtable(source, 'Delimiter', ',');
    fprintf('Loaded %d events from %s\n', height(T), source);
end

T.Properties.VariableNames = lower(strtrim(T.Properties.VariableNames));
require_columns(T, {'strike', 'dip', 'rake'}, source);

ranges = struct('strike', [0 360], 'dip', [0 90], 'rake', [-180 180]);
problems = {};
for name = {'strike', 'dip', 'rake'}
    col = name{1};
    values = T.(col);
    if ~isnumeric(values)
        problems{end + 1} = sprintf('''%s'' must be numeric', col); %#ok<AGROW>
        continue
    end
    lim = ranges.(col);
    bad = find(isnan(values) | values < lim(1) | values > lim(2));
    if ~isempty(bad)
        problems{end + 1} = sprintf('''%s'' must be between %g and %g; rows %s are not', ...
            col, lim(1), lim(2), mat2str(bad(:)' + 1)); %#ok<AGROW>
    end
end
if ismember('uncertainty_deg', T.Properties.VariableNames)
    bad = find(~(T.uncertainty_deg > 0));
    if ~isempty(bad)
        problems{end + 1} = sprintf('''uncertainty_deg'' must be positive; rows %s are not', ...
            mat2str(bad(:)' + 1));
    end
end
if ~isempty(problems)
    error('focal_stress:badCatalog', '%s: %s', source, strjoin(problems, '; '));
end

if ~ismember('event_id', T.Properties.VariableNames)
    ids = arrayfun(@(i) sprintf('EV%03d', i), (1:height(T))', 'UniformOutput', false);
    T = [table(ids, 'VariableNames', {'event_id'}), T];
else
    T.event_id = as_text_ids(T.event_id);
end
if numel(unique(T.event_id)) < height(T)
    error('focal_stress:duplicateIds', '%s: event_id values must be unique', source);
end
end
