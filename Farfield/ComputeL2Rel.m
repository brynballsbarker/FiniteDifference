
function errors = ComputeL2Rel(uscex, uscat)

global N m;
global r0;
global delta_theta;


maxerror = 0;
jmaxerror = 0;

the_sum_num = 0;
the_sum_den = 0;
for j=1:m+1
    add_num = r0*delta_theta*(abs(uscat(N,j))-abs(uscex(N,j)))^2;
    the_sum_num = the_sum_num + add_num;
    
    add_den = r0*delta_theta*(abs(uscex(N,j)))^2;
    the_sum_den = the_sum_den + add_den;
    
    error = add_num/add_den;
    if error > maxerror
        maxerror = error;
        jmaxerror = j;
    end
end

L2Relerror = (the_sum_num/the_sum_den)^(1/2);

errors = [L2Relerror, jmaxerror, maxerror];