function output = opA(dim,j,k)

% Output phase point operator A_jk for dimension d (prime).
% If k is not given, j must be a 2-component vector.

if nargin == 2
    j1 = mod(j(1),dim);
    j2 = mod(j(2),dim);
elseif nargin == 3
    j1 = mod(j,dim);
    j2 = mod(k,dim);
else
    disp('Error in opA.m: incorrect number of arguments.');
    return
end

if ~isprime(dim)
    disp('Error in opD.m: dimension d must be a prime.');
    return
end
   
A0 = zeros(dim);
for ii = 0:(dim-1)
for jj = 0:(dim-1)
    A0 = A0 + opD(dim,[ii,jj])/dim;
end
end

D = opD(dim,[j1,j2]);
output = D*A0*D';

end