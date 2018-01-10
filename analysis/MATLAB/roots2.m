% ROOTS2(A,B,C) returns solutions X to quadratic A*X^2 + B*X + C = 0.
function x = roots2 (a, b, c)
  q = -(b + sign(b)*sqrt(b^2 - 4*a*c))/2;
  x = [q/a c/q];