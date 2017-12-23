function output = Wigner2(rho,p,q)

% Returns qubit Wigner function for state rho at phase space
% coordinates p and q (takes value 0,1). Input p can be a vector.

if nargin == 1
    output = zeros(2);
    A11 = MUB2pro(0,1) + MUB2pro(1,1) + MUB2pro(2,1);
    A12 = MUB2pro(0,2) + MUB2pro(1,2) + MUB2pro(2,1);
    A21 = MUB2pro(0,1) + MUB2pro(1,2) + MUB2pro(2,2);
    A22 = MUB2pro(0,2) + MUB2pro(1,1) + MUB2pro(2,2);
    output(1,1) = (trace(A11*rho)-1)/2;
    output(1,2) = (trace(A12*rho)-1)/2;
    output(2,1) = (trace(A21*rho)-1)/2;
    output(2,2) = (trace(A22*rho)-1)/2;
    return
elseif nargin == 3
    i = p;
    j = q;
elseif nargin == 2
    if length(p) == 2
        i = p(1);
        j = p(2);
    else
        disp('Input p should be 2-component');
    end
end

if (i==0)&&(j==0)
    A = MUB2pro(0,1) + MUB2pro(1,1) + MUB2pro(2,1);
elseif (i==0)&&(j==1)
    A = MUB2pro(0,2) + MUB2pro(1,2) + MUB2pro(2,1);
elseif (i==1)&&(j==0)
    A = MUB2pro(0,1) + MUB2pro(1,2) + MUB2pro(2,2);
elseif (i==1)&&(j==1)
    A = MUB2pro(0,2) + MUB2pro(1,1) + MUB2pro(2,2);
else
    disp('Input p and q should take values 0 and 1.');
end

output = (trace(A*rho)-1)/2;

end
    