function match_vector = numeric_pattern(original, pattern)
     SIZE = length(original) - length(pattern);
     match = zeros(1, SIZE);
     for i=1:SIZE
         match(i) = all(original(i:i-1+length(pattern)) == pattern);
     end
     match_vector = find(match == 1);
 