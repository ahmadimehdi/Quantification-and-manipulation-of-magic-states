function output = MUB2pro(b,v)

% Generate MUB projector v = 1,2 in basis b = 0,1,2.

v = MUB2vec(b,v);
output = v*v';

end 