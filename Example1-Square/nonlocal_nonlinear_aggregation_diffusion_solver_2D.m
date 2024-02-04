function [rho, free_energy]=nonlocal_nonlinear_aggregation_diffusion_solver_2D(P_partition, T_partition, t_start, t_end, tau, num_control_volume, measure_control_volume, max_iteration_step, tolerence, time_interval)
%% PDE
% \rho_t = \nabla \cdot \Big[\rho \nabla \big(H'(\rho)+V+W*\rho \big) \Big]  in $\Omega \times [0, T]$, \\
% \rho \nabla \big(H'(\rho)+V+W*\rho \big) \cdot \bm{n} = 0 & on $\Gamma \times [0, T]$, \\
% \rho(\cdot,0)=\rho_0 & in $\Omega$.

%% annotation
%The problem domain is square.

%% basis_type_rho==P1
matrix_size=[num_control_volume num_control_volume];
vector_size=num_control_volume;
P_basis_rho=P_partition;
T_basis_rho=T_partition;

%% Initialize the iteration in time t = tau.
rho_old_time=get_initial_vector('function_rho_exact', P_partition, 0);

%% Free energy at T = 0
FE = function_free_energy(P_partition, num_control_volume, rho_old_time, measure_control_volume);
fprintf('The free energy at T=0: %d*****\n', FE);

%% Iteration in time.
nn = 2;
num_time = int16(( (t_end - t_start) / time_interval ) + 1);
free_energy = zeros(num_time, 1);
free_energy(1) = FE;        %%The free energy at initial time
N=(t_end-t_start)/tau;
for n=0:N-1           
    current_time=t_start+tau*(n+1);
%     fprintf('*********************Current_time: %d*********************\n',current_time);
    
    %% Assemble the load vector of rho. 
    B=assemble_vector_2D_time(rho_old_time, current_time, P_partition, vector_size, num_control_volume, measure_control_volume, tau);
    
    %% Iteration method
     %%Initialize the iteration
    rho_old_n=rho_old_time;    
    for l=1:max_iteration_step
%         fprintf('*********************Iteration %d *********************\n',l);
        
        %% Assemble the matrix A.   
        A=assemble_matrix_nonlinear_2D(P_partition, T_partition, tau, matrix_size, num_control_volume, measure_control_volume, rho_old_n, rho_old_time);
        rho_old=A\B;        
        err_iteration = norm(rho_old_n-rho_old);
%         fprintf('Iteration error %d\n',err_iteration);
%         rho_T_total_mass=rho_total_mass(rho_old, measure_control_volume, num_control_volume);
%         fprintf('******Current_time T=%d: the tatal mass %d at %d iteration*****\n',current_time, rho_T_total_mass, l);
%         if l == 1
%             rho_T0_total_mass=rho_total_mass(rho_old, measure_control_volume, num_control_volume);
%             fprintf('******Current_time T=%d: the tatal mass %d at %d iteration*****\n',current_time, rho_T0_total_mass, l);
%         end
        if (err_iteration)<tolerence    
%             fprintf('********Iteration %d step to convergence!*******\n',l);
%             rho_T_total_mass=rho_total_mass(rho_old, measure_control_volume, num_control_volume);
%             fprintf('******Current_time T=%d: the tatal mass %d at %d iteration*****\n',current_time, rho_T_total_mass, l);
            break;
        end
        if l == max_iteration_step
            fprintf('Not convergence!!!!!!!!!\n');
        end
        rho_old_n=rho_old;
    end
    rho_old_time=rho_old;  
    
    if rem(current_time,time_interval)==0
        %% Free energy
        FE = function_free_energy(P_partition, num_control_volume, rho_old_time, measure_control_volume);
        free_energy(nn) = FE;
        fprintf('Current_time T=%d: the free energy %d\n',current_time, FE);
        nn = nn + 1;
    
    end
end
rho=rho_old_time;

