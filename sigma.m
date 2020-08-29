%This function calculates the product of the raising operator of atom i
%with the lowering operator of atom j
%sigma^{i}_{eg}\times sigma^{j}_{ge}

function mat = sigma(i,j,N)
seg = sparse([0 1;0 0]);
sge = sparse([0 0;1 0]);
if i == j 
    seg = sparse([1 0;0 0]);
end
if i == 1
    mat = seg;
elseif j == 1
    mat = sge;
else
    mat = speye(2);
end

for it = 2:N
    if it == i
        matIt = seg;
    elseif it == j
        matIt = sge;
    else
        matIt = speye(2);
    end
    
    mat = kron(mat,matIt);
end

end

