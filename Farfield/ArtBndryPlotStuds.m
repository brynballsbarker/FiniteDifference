
function ArtBndryPlotStuds(NTAE,uscex,uscat,jmaxerr,maxerror,RelL2error)

% NTAE is the number of terms included in the Karp Expansion
% uscaex(i,j) is the exact solution matrix
% uscat(i,j) is the numerical solution matrix
% jmaxerr is the angular index where the maximum error occurs
% maxerror is the value of the maximunm error
% RelL2error is the L^2-Norm Rel Error made by approximating uscat by uscex

%Globals
global k rmax;
global theta;
global N m;


% %Comparison Plot
figure
plot(theta(1:m+1),abs(uscex(N,1:m+1)),'k',theta(1:m+1),abs(uscat(N,1:m+1)),'r*')
str={['Comparison of Scattered Fields: '...
    'Analytical vs Approximated at Artificial Bdry  ']...
     ,['   Terms Expansion = ' num2str(NTAE)],...
    [num2str(N) ' x '  num2str(m+1)   '    k= ' num2str(k) '      R = '...
    num2str(rmax)],['MaxError = ' num2str(maxerror), '   at \theta = ' ...
    num2str(theta(jmaxerr)), '      L^2-Norm Rel Error = ' num2str(RelL2error)]};
xlabel('theta');
ylabel('Scattered field');
legend('Analytical', 'Numerical');
title(str);
hold on;
plot(theta(jmaxerr),abs(uscat(N,jmaxerr)), 'ko');

end