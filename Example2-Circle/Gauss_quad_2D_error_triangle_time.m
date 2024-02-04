function r=Gauss_quad_2D_error_triangle_time(rho_local,accurate_function,t,vertices,basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,derivative_degree_x,derivative_degree_y, a)
%Use Gauss quadrature to numerically compute a norm error of FE solution on a local triangular element T.
%accurate_function: the accurate function in the error.
%When we take the L2 norm,accurate_function is the exact solution.
%When we take the H1 seminorm, accurate_function is the first derivative of the exact solution.
%vertices: the coordinates of the vertices of the triangular element T.
%rho_local: the values of the FE solution at the nodes of FE in the triangular element T.
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangle.
%derivative_degree_x:the derivative degree of the FE solution with respect to x.
%derivative_degree_y:the derivative degree of the FE solution with respect to y.
%Gpn: the Gauss point number.

Gpn=length(Gauss_coefficient_reference_triangle);

r=0;
[Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);

for i=1:Gpn
    r=r+Gauss_coefficient_local_triangle(i)*(feval(accurate_function,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),t)-a*FE_solution_triangle_index(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),rho_local,vertices,basis_index,basis_type,derivative_degree_x,derivative_degree_y))^2;
end