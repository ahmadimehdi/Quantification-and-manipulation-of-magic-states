% Generate random density matrices in dimension d = 2 and check if there
% exists conversion between 2 states that increases mana.

d = 2; % dimension
n = 8; % number of phase point operators A_j
m = 6; % number of MUB projectors P_k

% Pauli operators
I = eye(2);
S(:,:,1) = [0,1;1,0]; 
S(:,:,2) = [0,-1i;1i,0];
S(:,:,3) = [1,0;0,-1];

A = zeros(d,d,n); % define phase point operators
for jj = 1:4
    A(:,:,jj) = opA(d,mod(floor((jj-1)/2),2),mod(jj-1,2));
    A(:,:,jj+4) = transpose(A(:,:,jj));
end

P = zeros(d,d,m); % define MUB projectors
for kk = 1:m
    P(:,:,kk) = MUB2pro(mod(floor((kk-1)/2),3),mod(kk-1,2)+1);
end

% create random density matrices
N = 10; % re-run or increase this number if all trials fail
Rho = zeros(2,2,N);
Mana = zeros(1,N);
for ii = 1:N
    v = 2*rand(3,1)-1;
    v = v/norm(v);
    v = v*rand;
    Rho(:,:,ii) = B2D(v);
    Mana(ii) = sum(sum(abs(Wigner2(Rho(:,:,ii)))));
end

tic
flag = 0;
for i1 = 1:N
for i2 = (i1+1):N
    if Mana(i2) - Mana(i1) > 1e-7
        rho1 = Rho(:,:,i1);
        rho2 = Rho(:,:,i2);
        mana1 = Mana(i1);
        mana2 = Mana(i2);
    elseif Mana(i1) - Mana(i2) > 1e-7
        rho2 = Rho(:,:,i1);
        rho1 = Rho(:,:,i2);
        mana2 = Mana(i1);
        mana1 = Mana(i2);
    else
        break;
    end

    cvx_begin quiet
        variable J(4,4) complex semidefinite % completely positive 
        subject to
        trace(J) == 2
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
        for ii = 1:3
            Si = S(:,:,ii);
            trace(J*kron(I,Si)) == 0
            trace(J*(kron(Si,transpose(rho1))-trace(Si*rho2)*kron(I,I)/2)) == 0
          end
    cvx_end
    if strcmp(cvx_status,'Solved')
        disp('Conversion detected.');
        disp('The following density matrix')
        disp(rho1);
        disp('with mana')
        disp(mana1);
        disp('can be converted to:');
        disp(rho2);
        disp('with mana ')
        disp(mana2)
        disp('by the following Choi matrix')
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
if flag == 0 
    disp('No conversion found, please run again.')
end