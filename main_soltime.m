%% Robust MPC with Parametric Uncertainty: Main Code to get Time Values 
% Solves Nmax MPC problems with horizons 1 to Nmax from Ninit samples
% initial conditions

clear all
close all
clc
yalmip 'clear'

%% Matrix Uncertainty Parameters
epsA = 0.1;     
epsB = 0.1;      % These are infinity norm bounds on the matrix uncertainties. 

%% Load all required sets and parameters for MPC
[Anom,Bnom, delAv, delBv, K, A, B, X, U, Xlb, Xub, Ulb, Uub, nx, nu, wub, wlb, Q, R, N_max, N_thres] = sys_load(epsA, epsB);
W = Polyhedron('lb',wlb*ones(nx,1),'ub',wub*ones(nx,1));

%% Create sample initial conditions to solve from   
Ninit = 100; 
vec1 = linspace(Xlb(1), Xub(1),sqrt(Ninit)); 
vec2 = linspace(Xlb(2), Xub(2),sqrt(Ninit)); 
x_s = []; 
for i = 1:length(vec1)
    for j = 1:length(vec2)
        x_s = [x_s, [vec2(j); vec1(i)]]; 
    end
end 

%% Form the terminal set and Cost here 
[Xn, Pinf] = term_setRobPar(Anom, Bnom, delAv, delBv, K, X, U, W, Q, R, nx, nu); 

%% Get all possible vertices of matrix uncertainty
% Needed for constraint loop
for i = 1:size(delAv,2)/nx
    setdelA(:,:,i) = delAv(:,(i-1)*nx + 1: i*nx);  
end

for i = 1:size(delBv,2)/nu
    setdelB(:,:,i) = delBv(:,(i-1)*nu + 1: i*nu);  
end

%% Getting all the bounds and offline time values
bound_time = zeros(N_max,1); 
for Nhor = 1:N_max     
    Fx = blkdiag(kron(eye(Nhor-1), X.A), Xn.A); 
    boldAvbar = obtain_boldAvbar(Nhor, nx);        
    tic
    [t_w{Nhor}, t_1{Nhor}, t_2{Nhor}, t_3{Nhor}, t_delTaA{Nhor}, t_delTaB{Nhor}] = bounds(Fx, Anom, Bnom, Nhor, N_thres, boldAvbar, delAv, delBv, nx, nu);
    bound_time(Nhor) = toc; 
end

%% Getting Started with Nmax MPC Problems from each Sample State                                               
% Record all the time values for all horizon lengths and take an average
for i=1:Ninit
    x_init = x_s(:,i); 
    % solve with varying horizons here. Pick the best cost and go ahead 
    for Nhor = 1:N_max 
        [feas_flag(Nhor), cost_flag(Nhor), v_horN{Nhor}, sol_time(Nhor, i)] = FTOCP_Time(x_init, Nhor, Anom, Bnom, Xn, X, U, W, wub, nx, nu, Q, R, Pinf, ...
                                                                setdelA, setdelB, t_w{Nhor}, t_1{Nhor}, t_2{Nhor}, t_3{Nhor}, t_delTaA{Nhor}, t_delTaB{Nhor}); 
    end
    yalmip 'clear'         
end

%% Get the mean solver time 
online_time = mean(sol_time,2);                                                    % for each horizon

