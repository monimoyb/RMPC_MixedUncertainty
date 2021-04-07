%% Get Approximate ROA with Optimization
clear all
close all
clc
yalmip 'clear'

%% Matrix Uncertainty Parameters
vertflag = 0;    % Use 1 for full vertices required to represent a norm uncertainty. Use 0 for structured uncertainty. 
epsA = 0.1;     
epsB = 0.1;      % These are infinity norm bounds on the matrix uncertainties. 

%% Load all required sets and parameters for MPC
[Anom,Bnom, delAv, delBv, K, A, B, X, U, Xlb, Xub, Ulb, Uub, nx, nu, wub, wlb, Q, R, N_max, N_thres] = sys_load(vertflag, epsA, epsB);
W = Polyhedron('lb',wlb*ones(nx,1),'ub',wub*ones(nx,1));

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

%% Compute Bounds
if N_max == 1
   Fx = Xn.A; 
   t_w{1} = zeros(size(Fx,1),1);
   t_1{1} = zeros(size(Fx,1),1);
   t_2{1} = zeros(size(Fx,1),1); 
   t_3{1} = zeros(size(Fx,1),1);
else
%% Form the other system matrices and load all the bounds here 
  Fx = blkdiag(kron(eye(N_max-1), X.A), Xn.A); 
  boldAvbar = obtain_boldAvbar(N_max, nx);        
  [t_w{1}, t_1{1}, t_2{1}, t_3{1}] = bounds(Fx, Anom, Bnom, N_max, N_thres, boldAvbar, delAv, delBv, nx, nu);
end                                             

%% Compute the ROA i.e., Nmax-Step Rob. Controllable Set
%%% Add as many direction vectors as you wish
dVector{1} =  [1;1];
dVector{2} = [ 0;1];
dVector{3} = [ 1;0];
dVector{4} = [-1;1];
dVector{5} = [ 2;-6];
dVector{6} = [2;6];
dVector{7} = [-6;8.2];
dVector{8} = [8.1;-6.2];
dVector{9} = [8.0;-4.439];
%% search in the negative directions too
vSign{1}    =  1;
vSign{2}    = -1;

%% Main Loop 
x0feas = [];

for i = 1:length(dVector)
    for j = 1:2
        [x0feas_out, x0feasNormOut] = FTOCP(dVector{i}, vSign{j}, N_max, Anom, Bnom, Xn, X, U, W, wub, nx, nu, ...
                                                                                   setdelA, setdelB, t_w{1}, t_1{1}, t_2{1}, t_3{1});   
        if x0feasNormOut ~= -inf
            x0feas = [x0feas, x0feas_out];   % only if feasible 
        end
    end
end

%% Plot ROA i.e., Nmax-step rob. controllable set
feasPolMPC_app = Polyhedron(x0feas'); 
figure;  
plot(feasPolMPC_app,'color','b','alpha',1); hold on; 
set(gca, 'fontsize',25,'fontweight','bold');
xlim([-8.2,8.2]); ylim([-8.2,8.2]);
