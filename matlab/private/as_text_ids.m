function ids = as_text_ids(ids)
%AS_TEXT_IDS Return event IDs as a column cell array of char.
if isnumeric(ids)
    ids = arrayfun(@(x) sprintf('%g', x), ids, 'UniformOutput', false);
elseif isa(ids, 'string')
    ids = cellstr(ids);
end
ids = ids(:);
end
