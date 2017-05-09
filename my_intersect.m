function C = my_intersect(A,B)
% MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
% C = my_intersect(A,B)

% by Leon Peshkin pesha at ai.mit.edu, after Kevin Murphy's BNT
A = A(:)'; B = B(:)';
if isempty(A) | isempty(B)
  C = [];
  return
else
  bits = zeros(1, max(max(A),max(B)));  %bits = sparse(1, max(ma,mb));
  bits(A) = 1;
  C = B(logical(bits(B)));  
end