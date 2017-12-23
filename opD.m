function output = opD(dim,p)

% Create Weyl-Heisenberg displacement operator D.
% If only 2 inputs are given (d and p1), then p1 will 
% be understood as vector p1 that runs from 1 to d^2.

if nargin == 1                  % if only input p is given
    d = getd;    
elseif nargin == 2
    d = dim;
end

if length(p) == 1               % single index
        p2 = mod(p,d);
        p1 = mod(floor(p/d),d);
elseif length(p) == 2           % double index
        p2 = p(2);
        p1 = p(1);
end

if dim > 2
    tau = exp(pi*1i*(d+1)/d);               
    output = tau^(p1*p2)*(opX(d)^p1)*(opZ(d)^p2);
elseif dim == 2
    output = 1i^(p1*p2)*(opX(d)^p1)*(opZ(d)^p2);
end

end