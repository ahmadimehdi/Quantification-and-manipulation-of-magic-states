% Given qubits rho1 and rho2, determine if there exists
% a free channel converting rho1 to rho2.

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

% % create random points inside Bloch sphere
% v1 = rand(3,1);
% v1 = rand*v1/norm(v1);
% v2 = rand(3,1);
% v2 = rand*v2/norm(v2);
% 
% rho1 = B2D(v1);
% rho2 = B2D(v2);

% create random density matrices
N = 10; % number of vectors
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
    elseif Mana(i1) - Mana(i2) > 1e-7
        rho2 = Rho(:,:,i1);
        rho1 = Rho(:,:,i2);
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
toc


% cvx_begin quiet
%     variable J(4,4) complex semidefinite % completely positive 
%     subject to
%     trace(J) == 1
%     % conditions for free channels
%     for jj = 1:n
%     for kk = 1:m
%         Aj = A(:,:,jj);
%         Pk = P(:,:,kk);
%         PkT = transpose(Pk);
%         real(trace(J*kron(Aj,PkT))) >= 0
%     end
%     end
%     % conditions for trace preserving and state conversion
%     for ii = 1:3
%         Si = S(:,:,ii);
%         trace(J*kron(I,Si)) == 0
%         trace(J*(kron(Si,transpose(rho1))-trace(Si*rho2)*kron(I,I)/2)) == 0
%       end
% cvx_end
% if strcmp(cvx_status,'Solved')
%     disp('Conversion is possible with the following Choi matrix:');
%     disp(J);
% else
%     disp('Conversion is not possible.');
% end
