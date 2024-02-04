function measure_control_volume = function_measure_control_volume(P_partition, T_partition, num_control_volume)

%% Initialization
measure_control_volume = zeros(num_control_volume, 1);

%% 
for i = 1:num_control_volume
    [m,n] = find(T_partition(1:3, :) == i);    
    num_triangle = length(n);
    vertices_A =  P_partition(:,i); %%��ǰi������Ľڵ�����
    measure_control_volume_i = 0;   %%������i�����
    
    for j = 1:num_triangle
        a = T_partition(1:3, n(j));     %%ȷ��i��������Χ�ĵ�j�������νڵ���
        a(a == i) = [];                 %%�õ�i��������Χ�ĵ�j�������ε����������ڵ���        
        vertices_B =  P_partition(:,a(1)); %%��ǰi��������Χ�ĵ�j�������εĵ�һ���ڵ�����
        vertices_C =  P_partition(:,a(2)); %%��ǰi��������Χ�ĵ�j�������εĵڶ����ڵ�����
        vertices_AB_mid = (vertices_A + vertices_B)/2;  %%������AB�ߵ��е�����
        vertices_AC_mid = (vertices_A + vertices_C)/2;  %%������AC�ߵ��е�����
        vertices_tri_ij = circumcenter(vertices_A(1), vertices_A(2), vertices_B(1), vertices_B(2), vertices_C(1), vertices_C(2));      %%�����ε���������
        %% ���������μ���������
%             x = [vertices_A(1), vertices_B(1), vertices_C(1), vertices_tri_ij(1)];
%             y = [vertices_A(2), vertices_B(2), vertices_C(2), vertices_tri_ij(2)];
%             plot(x, y, 'r*');
%             axis([-1, 1, -1, 1]);
        %%
        edge_A_ABmid = sqrt((vertices_AB_mid(1) - vertices_A(1))^2 + (vertices_AB_mid(2) - vertices_A(2))^2);    %%AB���е㵽i������ľ���(����ֻ�����AB�߳��ȵ�һ��)
        edge_A_ACmid = sqrt((vertices_AC_mid(1) - vertices_A(1))^2 + (vertices_AC_mid(2) - vertices_A(2))^2);    %%AC���е㵽i������ľ���(����ֻ�����AC�߳��ȵ�һ��)
        edge_midAB_circ = sqrt((vertices_AB_mid(1) - vertices_tri_ij(1))^2 + (vertices_AB_mid(2) - vertices_tri_ij(2))^2);    %%AB���е㵽���ĵľ���
        edge_midAC_circ = sqrt((vertices_AC_mid(1) - vertices_tri_ij(1))^2 + (vertices_AC_mid(2) - vertices_tri_ij(2))^2);    %%AC���е㵽���ĵľ���
        measure_control_volume_i =  measure_control_volume_i + 1/2 * edge_A_ABmid * edge_midAB_circ + 1/2 * edge_A_ACmid * edge_midAC_circ;   %%������i�Ĳ��
    end
    measure_control_volume(i) = measure_control_volume_i;      %%i������Ĳ��    
end
