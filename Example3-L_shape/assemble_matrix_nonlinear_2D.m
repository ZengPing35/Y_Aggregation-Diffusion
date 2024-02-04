function A=assemble_matrix_nonlinear_2D(P_partition, T_partition, tau, matrix_size, num_control_volume, measure_control_volume, rho_old_n, rho_old_time)
%%�� Up-wind scheme ����װ $\rho_t = \nabla \cdot \Big[\rho \nabla \big(H'(\rho)+V+W*\rho \big) \Big]$
%% Initialization
A = sparse(matrix_size(1),matrix_size(2)); 
function_discrete_convolution_vector=zeros(num_control_volume,1);
xi_vector=zeros(num_control_volume,1);
grad_H=function_vector_grad_H(rho_old_n);

%% Function convolution; Function grad_H; Function xi
for n = 1:num_control_volume
    vertices=P_partition(:,n);              
    function_discrete_convolution_vector(n)=function_Discrete_convolution(P_partition, num_control_volume, vertices, measure_control_volume, rho_old_n, rho_old_time);
    xi_vector(n) = grad_H(n) + function_V(vertices(1,1), vertices(2,1)) + function_discrete_convolution_vector(n);
end

%% Up-wind scheme
for i = 1:num_control_volume
    A(i, i) = measure_control_volume(i) / tau;  %%rho_t��һ�������
    xi_i = xi_vector(i);                        
    vertices_A =  P_partition(:,i);             %%��ǰi������Ľڵ�����   
    [m,n] = find(T_partition == i);             %%n���������i��Χ�����εı��
    b = T_partition(:, n);                      %%ȷ��i��������Χ���������νڵ���
    c = b(:);               %%����������
    c(c == i) = [];         %%ȥ��i��
    d = unique(c);          %%ȥ���ظ�Ԫ��
    num = length(d);  
    for j = 1:num                   %%i��������Χ������ĸ���
        xi_j = xi_vector(d(j));     
        [q,t] = find(b == d(j));    %%�жϵ�j������������ڿ�����i��Χ���ļ�����������
        num_vol_in_tri = length(t); 
        vertices_B =  P_partition(:,d(j));  %%��j������������ڵ�k�������εĵڶ����ڵ�����
        vertices_AB_mid = (vertices_A + vertices_B)/2;  %%������AB�ߵ��е�����
        edge_AB = sqrt((vertices_B(1) - vertices_A(1))^2 + (vertices_B(2) - vertices_A(2))^2);    %%AB��ľ���
        edge_ij = 0;                        %%��ʼ��edge_ij
        for k = 1:num_vol_in_tri
            e = b(:, t(k));                 %%��j������������ڵ�k�������ε������νڵ���
            e(e == i) = [];                 %%��k�������γ�ȥi�ڵ� 
            e(e == d(j)) = [];              %%��ȥ����k�������εĽڵ�
            vertices_C =  P_partition(:,e); %%��j������������ڵ�k�������εĵ������ڵ�����                   
            vertices_tri_ij = circumcenter(vertices_A(1), vertices_A(2), vertices_B(1), vertices_B(2), vertices_C(1), vertices_C(2));      %%�����ε���������
            
           %% ���������μ���������
%             x = [vertices_A(1), vertices_B(1), vertices_C(1), vertices_tri_ij(1)];
%             y = [vertices_A(2), vertices_B(2), vertices_C(2), vertices_tri_ij(2)];
%             plot(x, y, 'r*');
%             axis([-1, 1, -1, 1]);
            
           %% edge:(AB�е㵽�����߶�)--->�����AB�е㵽�Ա��������������ĳ���֮��           
            edge_ij = edge_ij + sqrt((vertices_tri_ij(1) - vertices_AB_mid(1))^2 + (vertices_tri_ij(2) - vertices_AB_mid(2))^2); %%AB���е㵽���ĵľ���
        end
        upwind = (edge_ij / edge_AB) * (xi_i - xi_j);
        if upwind >= 0
            A(i, i) = A(i,i)+upwind;
        else
            A(i, d(j)) = upwind;
        end
    end
end