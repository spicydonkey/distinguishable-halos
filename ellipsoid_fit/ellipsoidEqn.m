function fellipsoid = ellipsoidEqn(v,x,y,z)
% F = ELLIPSOIDEQN(PV,X,Y,Z)
%
% Evaluates the fully parameterised ellipsoid function:
% F_ELLIP = Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J
%
% V:    10 parameters definining ellipsoid
% X,Y,Z:    cartesian grid points to evaluate function
% F:    ellipsoid function evaluated at the grid points
%
% NOTE: ellipsoid is the 2D root surface of the function
%

fellipsoid = v(1)*x.*x +   v(2)*y.*y +   v(3)*z.*z + ...
           2*v(4)*x.*y + 2*v(5)*x.*z + 2*v(6)*y.*z + ...
           2*v(7)*x    + 2*v(8)*y    + 2*v(9)*z    + v(10);
       
end