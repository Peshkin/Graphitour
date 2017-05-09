function b = my_unique(a) 
%UNIQUE Set unique.
%   UNIQUE(A) for the array A returns the same values as in A but
%   with no repetitions.  A will also be sorted.  A can be a cell
%   array of strings.

% by Leon Peshkin, inspired by BNT of Kevin Murphy

b = sort(a);
  % d indicates the location of matching entries
d = b((1:end-1)') == b((2:end)');
b(find(d)) = []; 