% Given qutrits rho1 and rho2, determine if there exists
% a free channel converting rho1 to rho2.

d = 3; % dimension
n = 81; % number of phase point operators A_j
m = 12; % number of MUB projectors P_k

I = eye(d);
% basis for traceless Hermitian
S(:,:,1) = [1,0,0;0,0,0;0,0,-1]; 
S(:,:,2) = [0,0,0;0,1,0;0,0,-1];
S(:,:,3) = [0,1,0;1,0,0;0,0,0];
S(:,:,4) = [0,0,1;0,0,0;1,0,0];
S(:,:,5) = [0,0,0;0,0,1;0,1,0];
S(:,:,6) = [0,1i,0;-1i,0,0;0,0,0];
S(:,:,7) = [0,0,1i;0,0,0;-1i,0,0];
S(:,:,8) = [0,0,0;0,0,1i;0,-1i,0];

% define phase point operators A_j
A = zeros(d,d,d^(d+1));
jj = 0;
for v1 = 1:d
for v2 = 1:d
for v3 = 1:d
for v4 = 1:d
    jj = jj + 1;
    Atemp = zeros(d,d);
    v = [v1, v2, v3, v4];
    for bb = 1:(d+1)
        Atemp = Atemp + MUB3pro(bb,v(bb));
    end
    A(:,:,jj) = Atemp - eye(d);
end
end
end
end

% MUB projectors P_k
P = zeros(d,d,m);
kk = 0;
for bb = 1:(d+1)
for vv = 1:d    
    kk = kk + 1;
    P(:,:,kk) = MUB3pro(bb,vv);
end
end

% create random density matrices
% N = 10;
% Rho = zeros(d,d,N);
% Mana = zeros(1,N);
% for ii = 1:N
%     Rho(:,:,ii) = randRho(d);
%     Mana(ii) = sum(sum(abs(Wigner(Rho(:,:,ii)))));
% end

% test for T and H state
N = 2;
v1 = [0;1;-1]/sqrt(2);
rho1 = v1*v1';
v2 = [-1;2;-1]/sqrt(6);
rho2 = v2*v2';

flag = 0;
for i1 = 1:N
for i2 = (i1+1):N
    if Mana(i2) - Mana(i1) > 1e-7
        rho1 = Rho(:,:,i1);
        rho2 = Rho(:,:,i2);
    elseif Mana(i1) - Mana(i2) > 1e-7
        rho2 = Rho(:,:,i1);
        rho1 = Rho(:,:,i2);
    else
        break;
    end
    
    tic
    cvx_begin quiet
        variable J(d^2,d^2) complex semidefinite % completely positive 
        subject to
        trace(J) == d
        % conditions for free channels
        for jj = 1:n
        for kk = 1:m
            Aj = A(:,:,jj);
            Pk = P(:,:,kk);
            PkT = transpose(Pk);
            real(trace(J*kron(Aj,PkT))) >= 0
        end
        end
        % conditions for trace preserving and state conversion
        for ii = 1:(d^2-1)
            Si = S(:,:,ii);
            trace(J*kron(I,Si)) == 0
            trace(J*(kron(Si,transpose(rho1))-trace(Si*rho2)*kron(I,I)/d)) == 0
          end
    cvx_end
    toc
    if strcmp(cvx_status,'Solved')
        [i1,i2]
        disp('Conversion is possible with the following Choi matrix:');
        disp(J);
        flag = 1;
        break;
    else
        disp('Conversion is not possible.');
    end
end
if flag == 1
    break
end
end