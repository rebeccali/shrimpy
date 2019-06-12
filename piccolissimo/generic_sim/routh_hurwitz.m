function result_array = routh_hurwitz(coeff)
n = length(coeff);
rows = n;
for cols = 1:ceil(rows/2)
    result_array(n-rows+1,cols) = coeff(2*(cols-1)+1);
    if strcmp(class(coeff),'sym')
        result_array(n-rows+1,cols) = simplify(result_array(n-rows+1,cols));
    end
end
rows = n-1;
for cols = 1:ceil(rows/2)
    result_array(n-rows+1,cols) = coeff(2*(cols-1)+2);
    if strcmp(class(coeff),'sym')
        result_array(n-rows+1,cols) = simplify(result_array(n-rows+1,cols));
    end
end
for rows = n-2:-1:1
    for cols = 1:ceil(rows/2)
        result_array(n-rows+1,cols) = -1/result_array(n-rows,1)*det([result_array(n-rows-1,1), result_array(n-rows,1);result_array(n-rows-1,cols+1), result_array(n-rows,cols+1)]);
        if strcmp(class(coeff),'sym')
            result_array(n-rows+1,cols) = simplify(result_array(n-rows+1,cols));
        end
    end
end
if strcmp(class(coeff),'sym')
    result_array = simplify(result_array,10000);
end