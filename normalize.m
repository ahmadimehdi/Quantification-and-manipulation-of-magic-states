function output = normalize(A);

% Normalizes a vector A. For matrix this will normalize
% in the sense of matrix 2-norm.

output = A/norm(A);

end