function output = MUBp(p)

% this construction is only for prime dimensions

output = zeros(p,p*(p+1));
output(:,1:p) = eye(p);
w = exp(2*pi*1i/p);
tau = w^((p+1)*inverse(2,p));

for kk = 1:(p-1)
for rr = 0:(p-1)
for ss = 0:(p-1)
   output(rr+1,kk*p+ss+1) = tau^(mod(inverse(kk,p)*(ss-rr)^2,p))/sqrt(p);
end
end
end

% Fourier basis
for rr = 0:(p-1)
for ss = 0:(p-1)
   output(rr+1,p^2+ss+1) = tau^(-mod(inverse(kk,p)*2*rr*ss,p))/sqrt(p);
end
end

end