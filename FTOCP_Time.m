%% Finite time optimal control problem: This one used to get the time values


function [feas_flag, cost_flag, v_hor, sol_time] = FTOCP_Time(x_0, N, Anom, Bnom, Xn, X, U, W, wub, nx, nu, Q, R, Pinf, setdelA, setdelB, t_w, t_1, t_2, t_3)

options = sdpsettings('solver','gurobi','verbose',0);

boldA1bar = zeros(nx*N, nx*N);
for j = 1:N
    tmpMat = [];
    for k = 1:j
        tmpMat = [tmpMat, Anom^(j-k)];
    end
    boldA1bar(nx*(j-1)+1: nx*j, 1:nx*j) = tmpMat;    
end
%% Need to vary N in a loop and solve the optimization problems. Have to pick the best cost
    boldAbar = kron(eye(N),Anom);
    boldBbar = kron(eye(N),Bnom); 
    %% Forming other matrices appearing in the optimization problem 
    Fx = blkdiag(kron(eye(N-1), X.A), Xn.A); 
    fx = [kron(ones(N-1,1),X.b); Xn.b]; 
    boldHw = kron(eye(N), W.A);
    boldhw = kron(ones(N,1),W.b); 
    boldHu = kron(eye(N), U.A);
    boldhu = kron(ones(N,1), U.b); 
%%%

    constraints = []; 
   %% Creating Open Loop Optimization Variables for MPC    
    M = sdpvar(nu*N,nx*N);       
        for j=1:nu*N
            for k = 2*j-1:nx*N
                    M(j,k) =0;
            end
        end
        
   v = sdpvar(nu*N,1);                                                     % nominal inputs
   xbf = sdpvar(nx*(N+1),1);                                               % nominal states 
   
   Lambda = sdpvar(size(Fx,1), size(boldHw,1), 'full');    
   gamma = sdpvar(size(boldHw,1), size(boldHu,1), 'full'); 
 
   %% Solving the Optimization Problems with varying horizon 
   xbf(1:nx,1) = x_0;        
   for k=1:N
        xbf(k*nx+1:(k+1)*nx, 1)  = Anom* xbf((k-1)*nx+1:k*nx, 1)  + Bnom*v(1+(k-1)*nu:k*nu,1);
   end
   
    cost_state =   xbf'*blkdiag(kron(eye(N),Q), Pinf)*xbf;
   
   %% Cleanly separate the cases here
   if N ==1 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% enumerate the set of vertices here 
    for ii = 1:size(setdelA,3)
        bolddelA = kron(eye(N),setdelA(:,:,ii));
        for jj = 1:size(setdelB,3) 
                bolddelB = kron(eye(N),setdelB(:,:,jj));
                constraints = [constraints; Fx*((boldAbar+bolddelA)*x_0 + (boldBbar+bolddelB)*v) + Lambda*boldhw...
                                                                                      <= fx];
        end
     end
   
    constraints = [constraints; Lambda*boldHw == Fx];
  
   else
     
 %%%% Now in this case I do not want to dualize Dela and Delb!!              
     constraints = [constraints; Fx*boldAbar*xbf(1:end-nx,1) + Fx*boldBbar*v...
                                     + (t_1)*norm(xbf(1:end-nx,1),inf) + (t_2)*(norm(M, inf)*wub+norm(v, inf)) + (t_3)*norm(M,inf)*wub...
                                     + (t_w)*wub + Lambda*boldhw <= fx];                                                                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     constraints = [constraints, Fx*boldBbar*M + Fx*(boldA1bar*boldBbar - boldBbar)*M + Fx*eye(size(boldA1bar,1)) == Lambda*boldHw]; 
      
   end
  
  %%%%%%%%%% These ones remain unaltered for all horizon lengths %%%%%%%%%%
  
   constraints = [constraints; Lambda>=0; gamma >=0];

   constraints = [constraints; gamma'*boldhw <= boldhu - boldHu*v];
      
   constraints = [constraints; (boldHu*M)' == boldHw'*gamma]; 
   
   obj_ol = v'*kron(R,N)*v + cost_state;   
   
   diagn=solvesdp(constraints, obj_ol, options);
   
   feas_flag = diagn.problem; 
   sol_time = diagn.solvertime;
   
   if feas_flag ~=0
       cost_flag = inf;                       % store high cost if any issue
       v_hor = zeros(nu*N,1);                 % just store some dummy 0's if infeasible anyway 
   else
       cost_flag = double(obj_ol);            % store right cost if feasible  
       v_hor = double(v); 
   end
   
  
end
