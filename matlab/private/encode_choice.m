function key = encode_choice(choice)
%ENCODE_CHOICE Plane choice [1 2 2 1] -> '1|2|2|1' (text, so no tool reads it as a number).
key = strjoin(arrayfun(@(c) sprintf('%d', c), choice(:)', 'UniformOutput', false), '|');
end
