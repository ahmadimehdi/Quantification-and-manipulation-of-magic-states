% generate random density matrices and
% calculate their "mana"
N = 1000; % number of vectors
Rho = zeros(2,2,N);
listM = zeros(1,N);
for ii = 1:N
    v = 2*rand(3,1)-1;
    v = v/norm(v);
    v = v*rand;
    Rho(:,:,ii) = B2D(v);
    listM(ii) = sum(sum(abs(Wigner2(Rho(:,:,ii)))));
end

