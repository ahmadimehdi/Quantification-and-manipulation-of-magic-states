d = 3;
A = zeros(d,d,d^(d+1));
count = 0;
for v1 = 1:d
for v2 = 1:d
for v3 = 1:d
for v4 = 1:d
    count = count + 1;
    Atemp = zeros(d,d);
    v = [v1, v2, v3, v4];
    for ii = 1:(d+1)
        Atemp = Atemp + MUB3pro(ii,v(ii));
    end
    A(:,:,count) = Atemp - eye(d);
end
end
end
end

Gram = zeros(81,81);
for ii = 1:81
for jj = 1:81
    Gram(ii,jj) = trace(A(:,:,ii)*A(:,:,jj));
end
end