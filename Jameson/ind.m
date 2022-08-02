function i = ind(j,n)
k = @(j) (mod(j-1,n) + 1);
if j - 2 >= 1 && j+2 <= n
    i = j;
else  
    i = k(j);
end
    