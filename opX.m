function output = opX(d);

% Create a d x d matrix for the shift operator X

output = zeros(d,d);
output(1,d) = 1;

for ii = 2:d
    output(ii, ii-1) = 1;
end