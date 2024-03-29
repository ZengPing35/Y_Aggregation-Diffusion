function main
%The problem domain is L_shape.

%% Parameters
format short e
max_iteration_step=60;
tolerence=10^(-6);
t_start=0;
t_end=1.1;
time_interval = 0.05;        %%Output figure-interval

%% Saving results
diary('data.txt');

%% fix tau = 1/40000
tau=1/40000;
fprintf('Fix (tau=%d)\n',tau);
        
%% numerical solution (h = 1/24)

% Mesh information for partition
[P_partition,T_partition,Edge] = getMesh_Gmsh('L_shape_24.msh');
T_partition = T_partition(1:3, :);
num_nodes = size(P_partition,2);            %%The nodes number of the triangulation mesh
num_control_volume = num_nodes;             %%The number of control volumes

% Measure of control volume
measure_control_volume = function_measure_control_volume(P_partition, T_partition, num_control_volume);

% numerical solution                                 
[rho, free_energy, Total_mass]=nonlocal_nonlinear_aggregation_diffusion_solver_2D(P_partition, T_partition, t_start, t_end, tau, num_control_volume, measure_control_volume, max_iteration_step, tolerence, time_interval);

% figure of free energy
figure;
free_energy_x = 0:time_interval:t_end;
plot(free_energy_x', free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_24.fig')
close;     

% figure of total mass
filename = 'Mass_conservation_24.mat';
save(filename, 'Total_mass');

figure;
Total_mass_x = 0:time_interval:t_end;
plot(Total_mass_x', Total_mass, 'r');
xlabel('time');
ylabel('The total mass of \rho^n_h');
title('Mass conservation');
savefig('Mass_conservation_24.fig')
close;    
