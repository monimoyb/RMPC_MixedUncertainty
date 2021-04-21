%% Finite time optimal control problem solved to get the Approximate ROA
% Monimoy Bujarbaruah

function [x0feas_out, x0feasNormOut] = FTOCP(dVector, vSign, N, Anom, Bnom, Xn, X, U, W, wub, nx, nu, setdelA, setdelB, t_w, t_1, t_2, t_3, t_delTaA, t_delTaB)

    options = sdpsettings('solver','gurobi','verbose',0);
    boldA1bar = zeros(nx*N, nx*N);

    for j = 1:N
        tmpMat = [];
        for k = 1:j
            tmpMat = [tmpMat, Anom^(j-k)];
        end
        boldA1bar(nx*(j-1)+1: nx*j, 1:nx*j) = tmpMat;    
    end

    %% Forming matrices appearing in the optimization problem 
    boldAbar = kron(eye(N),Anom);
    boldBbar = kron(eye(N),Bnom); 
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

    gamma = sdpvar(size(boldHw,1), size(boldHu,1), 'full'); 

    x0feas = sdpvar(nx,1);

    %% Solving the Optimization Problems with varying horizon 
    xbf(1:nx,1) = x0feas;        
    for k=1:N
        xbf(k*nx+1:(k+1)*nx, 1)  = Anom* xbf((k-1)*nx+1:k*nx, 1)  + Bnom*v(1+(k-1)*nu:k*nu,1);
    end

   %% Cleanly separate the cases here. 
    if N ==1 
        Lambda = sdpvar(size(Fx,1), size(boldHw,1), 'full'); 
        constraints = [constraints; Lambda>=0];
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% enumerate the set of vertices here 
        for ii = 1:size(setdelA,3)
            bolddelA = kron(eye(N),setdelA(:,:,ii));
            for jj = 1:size(setdelB,3) 
                bolddelB = kron(eye(N),setdelB(:,:,jj));
                constraints = [constraints; Fx*((boldAbar+bolddelA)*x0feas + (boldBbar+bolddelB)*v) + Lambda*boldhw...
                                                                                    <= fx + soft_flg*epsilon];
            end
        end

        constraints = [constraints; Lambda*boldHw == Fx];

    else 
        Lambda = sdpvar(size(Fx,1), size(boldHw,1), 'full');                                % dual variables
        constraints = [constraints, Lambda >=0]; 

        %%% Verex enumeration step here for two terms linear in model mismatches
        for ii = 1:size(setdelA,3)
            bolddelA = kron(eye(N),setdelA(:,:,ii));
                for jj = 1:size(setdelB,3) 
                    bolddelB = kron(eye(N),setdelB(:,:,jj));
                    constraints = [constraints; Fx*boldAbar*xbf(1:end-nx,1) + Fx*boldBbar*v + Fx*boldA1bar*bolddelA*xbf(1:end-nx,1) + Fx*boldA1bar*bolddelB*v...
                                     + (t_1)*norm(xbf(1:end-nx,1),inf) + (t_2+t_delTaB)*(norm(M, inf)*wub) + t_2*norm(v, inf) + (t_3)*norm(M,inf)*wub...
                                     + (t_w)*wub + Lambda*boldhw <= fx];
                end
        end
                                                        
        constraints = [constraints, Fx*boldBbar*M + Fx*(boldA1bar*boldBbar - boldBbar)*M + Fx*eye(size(boldA1bar,1)) == Lambda*boldHw];

    end

    %%%%%%%%%% These ones remain unaltered for all horizon lengths %%%%%%%%%%

    constraints = [constraints; gamma >=0];
    constraints = [constraints; gamma'*boldhw <= boldhu - boldHu*v];
    constraints = [constraints; (boldHu*M)' == boldHw'*gamma];   
    constraints = [constraints; dVector(2)*x0feas(1) == dVector(1)*x0feas(2)]; 
    constraints = [constraints; X.A*x0feas <= X.b];

    cost = vSign*dVector'*x0feas;  
    diagn=solvesdp(constraints, cost, options);
    feas_flag = diagn.problem; 

    if feas_flag ~=0
        cost_flag = inf;                     % store high cost if any issue
        x0feas_out = double(x0feas);
        x0feasNormOut = -inf;                % infeasibility flag
    else
        cost_flag = double(cost);            % store right cost if feasible  
        v_hor = double(v); 
        x0feas_out = double(x0feas);
        x0feasNormOut = norm(x0feas_out,2);
    end

end
