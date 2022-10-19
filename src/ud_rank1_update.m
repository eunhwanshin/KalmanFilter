function [U, D] = ud_rank1_update(U, D, c, a)
% The Agee-Turner UD update:
%   \bar{U} \bar{D} \bar{U}' = U D U' + c a a'
% Ref.: Bierman, G. J. (1977). Factorization Methods for Discrete Sequential
% Estimation. Academic Press, N.Y., p. 55.
%
n = size(U, 1);

for j = n:-1:2,
  s = a(j);
  Dj = D(j) + c * s * s;
  b = c / Dj;
  beta = s * b;
  c = b * D(j);
  D(j) = Dj;

  for i=1:j-1,
    a(i) = a(i) - s * U(i, j);
    U(i,j) = U(i,j) + beta * a(i);
  end

end

D(1) = D(1) + c * a(1) * a(1);


