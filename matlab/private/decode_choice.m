function choice = decode_choice(key)
%DECODE_CHOICE '1|2|2|1' -> [1; 2; 2; 1]
choice = str2double(strsplit(char(key), '|'))';
end
