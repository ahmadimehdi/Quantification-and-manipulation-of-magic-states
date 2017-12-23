function output = opZ(d);

% Create a d x d matrix for the phase operator Z

w = exp(2*pi*i/d);                      % root of unity
output = zeros(d,d);

for ii = 1:d
    output(ii,ii) = w^(ii-1);
end