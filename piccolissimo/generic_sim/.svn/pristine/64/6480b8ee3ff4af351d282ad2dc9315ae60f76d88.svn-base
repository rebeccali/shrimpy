function outMat = makeMat(vars)

%vars = [2 4 4]; yeilds a matrix with the first variable having 2 levels,
%and the second two with 4 variables

outMat = ones(prod(vars),size(vars,2));
orders(1) = 1;
for h = 2:size(vars,2)
    orders(h) = vars(h-1)*orders(h-1);
end

for i = 1:prod(vars)
    tempNum = i;
    for j = size(vars,2):-1:1
        for k = 1:vars(j)
            if tempNum > orders(j)
                outMat(i,j) = outMat(i,j) + 1;
                tempNum = tempNum - orders(j);
            end
        end
    end
end
    