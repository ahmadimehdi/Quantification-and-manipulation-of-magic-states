% calculate monotone for given sigma and t

d = 2; % dimension
n = 8; % number of phase point operators A_j
m = 6; % number of MUB projectors P_k

t = 0.5;
I = eye(2);
sigma = B2D(normalize([1,1,1]));
% sigma = [1,0;0,0]
rho0 = I/2;
% rho0 = B2D(normalize([1,1,1])/sqrt(3))
rho1 = B2D(normalize([1,1,1]));
alpha = 0.5;

value = [];
for alpha = 0:0.05:1
rho = (1-alpha)*rho0 + alpha*rho1;

P = zeros(d,d,d+1,d); % define MUB projectors
for bb = 1:(d+1)
for vv = 1:d    
    P(:,:,bb,vv) = MUB2pro(bb-1,vv);
end; end

PP = zeros(d^2,d^2,d,d,d,d+1,d);
for v1 = 1:d
for v2 = 1:d
for v3 = 1:d
for cc = 1:(d+1)
for uu = 1:d
    PP(:,:,v1,v2,v3,cc,uu) = kron(P(:,:,1,v1),P(:,:,cc,uu)) + ...
                             kron(P(:,:,2,v2),P(:,:,cc,uu)) + ...
                             kron(P(:,:,3,v3),P(:,:,cc,uu));
end; end; end; end; end

PP = reshape(PP,[4,4,48]);

cvx_begin quiet
    variable p(48) nonnegative
    variable S(2,2) hermitian semidefinite
    minimize(trace(S))
    subject to
    Omega = kron(sigma,transpose(rho));
    for ii = 1:48
        Omega = Omega + t*p(ii)*PP(:,:,ii);
    end
    kron(I,S) - Omega == hermitian_semidefinite(4)
    sum(p) == 1
  
cvx_end

value = [value; cvx_optval]
end

plot(0:0.05:1, value)
