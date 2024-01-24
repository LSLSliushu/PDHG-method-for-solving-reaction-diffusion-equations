

function [y] = CH_mollifier_function(x, epsilon)

   y = 2*exp(-epsilon*epsilon ./ (x.*x)).*(x < 0);

end