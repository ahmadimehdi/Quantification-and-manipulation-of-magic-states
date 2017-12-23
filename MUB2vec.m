function output = MUB2vec(b,v)

% Generate MUB vector v (1,2) in basis b (0,1,2). Only for qubit.

if (b==0)&&(v==1)
    output = [1;0];
elseif (b==0)&&(v==2)
    output = [0;1];
elseif (b==1)&&(v==1)    
    output = [1;1i]/sqrt(2);
elseif (b==1)&&(v==2)    
    output = [1;-1i]/sqrt(2);
elseif (b==2)&&(v==1)    
    output = [1;1]/sqrt(2);
elseif (b==2)&&(v==2)    
    output = [1;-1]/sqrt(2);
else
    disp('Input b should be 0,1,2, and v should be 1,2.');
end

end
    