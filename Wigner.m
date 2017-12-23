function output = Wigner(v)

% Wigner distribution of state v, or density matrix rho
% Note: only for odd dimensions.

[a,b] = size(v);
if a == b
    d = a;
    rho = v;
elseif a == 1
    d = b;
    rho = v'*v;
elseif b == 1
    d = a;
    rho = v*v';
else
    disp('Error: input has to be either a vector, or a square matrix.');
    return
end

A0 = circshift(fliplr(eye(d)),[0,1]);
output = zeros(d);
for ii = 0:(d-1)
for jj = 0:(d-1)
    Aij = opD(d,[ii,jj])*A0*opD(d,[ii,jj])';
    x = trace(rho*Aij)/d;
    if abs(imag(x)) > 1e-10
        disp('Error: Wigner function is not real-valued!');
    end
    output(ii+1,jj+1) = real(x);    
end
end

end
    