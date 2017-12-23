function [answer Choi]=convertible(rho1,rho2)
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
      answer = 1;
      Choi=J;
    else
        answer = 0;
        Choi=zeros(4);
    end
end