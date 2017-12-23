function output = MUB3vec(b,v)

% Generate MUB vector v (1,2,3) in basis b (1,2,3,4). Only for qubit.
w=exp(2*1i*pi/3);
if (b==1)&&(v==1)
    output = [1;0;0];
elseif (b==1)&&(v==2)
    output = [0;1;0];
elseif (b==1)&&(v==3)    
    output = [0;0;1];
elseif (b==2)&&(v==1)    
    output = [1;1;1]/sqrt(3); 
elseif (b==2)&&(v==2)    
    output = [1;w;conj(w)]/sqrt(3);
elseif (b==2)&&(v==3)    
    output = [1;conj(w);w]/sqrt(3);
elseif (b==3)&&(v==1)    
    output = [w;1;1]/sqrt(3); 
elseif (b==3)&&(v==2)    
    output = [1;w;1]/sqrt(3);
elseif (b==3)&&(v==3)    
    output = [1;1;w]/sqrt(3);
elseif (b==4)&&(v==1)    
    output = [conj(w);1;1]/sqrt(3); 
elseif (b==4)&&(v==2)    
    output = [1;conj(w);1]/sqrt(3);
elseif (b==4)&&(v==3)    
    output = [1;1;conj(w)]/sqrt(3);
else
    disp('Input b should be 1,2,3,4 and v should be 1,2,3.');
end

end
    