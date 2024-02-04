function result=FVM_solution_semiH1error(rho_exact, rho, P_partition, T_partition, num_control_volume)
%Numerically compute a norm error of FE solution on the whole Circle domain.
%rho: the values of the FE solution at all nodes of FE in the whole domain. These values are in 1D index of nodes of FE.
%accurate_function: the accurate function in the error.

%% Parameters
result1=0;%numerator of a fraction
result2=0;%the denominator of a fraction

%% Go through all elements and accumulate the error on them.
for i = 1:num_control_volume
    vertices_A =  P_partition(:,i);             %%��ǰi������Ľڵ����� 
    rho_i = rho(i);
    rho_exact_i = rho_exact(i);
    [m,n] = find(T_partition == i);             %%n���������i��Χ�����εı��
    b = T_partition(:, n);                      %%ȷ��i��������Χ���������νڵ���
    c = b(:);               %%����������
    c(c == i) = [];         %%ȥ��i��
    d = unique(c);          %%ȥ���ظ�Ԫ��
    num = length(d);  
    for j = 1:num                   %%i��������Χ������ĸ���
        [q,t] = find(b == d(j));    %%�жϵ�j������������ڿ�����i��Χ���ļ�����������
        num_vol_in_tri = length(t); 
        vertices_B =  P_partition(:,d(j));  %%��j������������ڵ�k�������εĵڶ����ڵ�����
        rho_j = rho(d(j));
        rho_exact_j = rho_exact(d(j));
        vertices_AB_mid = (vertices_A + vertices_B)/2;  %%������AB�ߵ��е�����
        edge_AB = sqrt((vertices_B(1) - vertices_A(1))^2 + (vertices_B(2) - vertices_A(2))^2);    %%AB��ľ���
        edge_ij = 0;                        %%��ʼ��edge_ij
        for k = 1:num_vol_in_tri
            e = b(:, t(k));                 %%��j������������ڵ�k�������ε������νڵ���
            e(e == i) = [];                 %%��k�������γ�ȥi�ڵ� 
            e(e == d(j)) = [];              %%��ȥ����k�������εĽڵ�
            vertices_C =  P_partition(:,e); %%��j������������ڵ�k�������εĵ������ڵ�����                   
            vertices_tri_ij = circumcenter(vertices_A(1), vertices_A(2), vertices_B(1), vertices_B(2), vertices_C(1), vertices_C(2));      %%�����ε���������
            
           %% edge:(AB�е㵽�����߶�)--->�����AB�е㵽�Ա��������������ĳ���֮��           
            edge_ij = edge_ij + sqrt((vertices_tri_ij(1) - vertices_AB_mid(1))^2 + (vertices_tri_ij(2) - vertices_AB_mid(2))^2); %%AB���е㵽���ĵľ���
        end
        rho_exact_ij = rho_exact_i - rho_exact_j;
        rho_numerical_ij = rho_i - rho_j;
        result1 = result1 + (abs(rho_exact_ij - rho_numerical_ij))^2 * edge_ij / edge_AB;
        result2 = result2 + (abs(rho_exact_ij))^2 * edge_ij / edge_AB;
    end
end

%% semi-norm
result1=sqrt(result1);
result2=sqrt(result2);
result=result1/result2;