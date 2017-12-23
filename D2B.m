function output = D2B(psi)

% Return the Bloch sphere coordinates of qubit state psi.
% Input psi can be a pure state vector or a density matrix.

[a,b] = size(psi);

if (a==1)&&(b==2)
    rho = psi'*psi;
elseif (a==2)&&(b==1)
    rho = psi*psi';
elseif (a==2)&&(b==2)
    rho = psi;
else
    disp('Error in D2B.m: psi must be a 2-component vector or 2x2 matrix.');
    return
end
   
u = 2*real(rho(1,2));
v = 2*imag(rho(2,1));
w = rho(1,1)-rho(2,2);

output = real([u,v,w]);
    
end

