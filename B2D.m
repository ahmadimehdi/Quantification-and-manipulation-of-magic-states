function output = B2D(v)

% Return the density matrix given Bloch vector v.
% Input v can be a 3 or 4-component vector. If v has
% 3 components, the coeff of I is assume to be 1/2.

n = length(v);
v = reshape(v,n,1);

if n==3
    v = [1; v];
elseif n ~= 4
    disp('Error in B2D.m: Input vector should have 3 or 4 components.');
    return
end

I = eye(2);
X = [0,1;1,0];
Y = [0,-1i;1i,0];
Z = [1,0;0,-1];
output = (v(1)*I + v(2)*X + v(3)*Y + v(4)*Z)/2;

end

