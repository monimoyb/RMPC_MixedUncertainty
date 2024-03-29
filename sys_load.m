%% Function defining System and Constraint Matrices 
%  Monimoy Bujarbaruah

function [Anom,Bnom, delAv, delBv, K, A, B, X, U, Xlb, Xub, Ulb, Uub, nx, nu, wub,wlb, Q, R, N_max, N_thres] = sys_load(epsA, epsB)

    %% Considering two states and one scalar input 

    Anom = [1, 0.15; 0.1, 1];
    Bnom = [0.1; 1.1];                    
    nx = size(Anom,2); nu = size(Bnom,2); 
    
    delAv = [ [0, epsA; epsA, 0], [0, epsA; -epsA, 0], [0, -epsA; epsA, 0], [0, -epsA; -epsA, 0]];     
    delBv = [[0; -epsB], [0; epsB], [epsB; 0], [-epsB; 0] ];

    %%%%%%%%% choose the feedback gain K here%%%
    K = place(Anom, Bnom, [0.745; 0.75]);             
    
    %% Set the true A and B matrices (satisfy the above bounds)
    A = [1, 0.05; 
           0, 1.0 ]; 
    B = [0; 
         1.1];                                                                                                                       

    %% Weights
    Q =  10*eye(nx);
    R =   2*eye(nu);

    %% Horizon and cut-off for bounds
    N_max = 3;             % max horizon N in the paper. Go up to 4 for good results 
    N_thres = 3;           % beyond this value do the binomial tail. 

    %% Considering constraints of the form -a<=x(i)<=a and -ulb<=u<=uub
    Xlb = -[8; 8];
    Xub = -Xlb; 
    Ulb =  -4; 
    Uub = -Ulb; 

    X = Polyhedron('lb',Xlb,'ub',Xub);
    U = Polyhedron('lb',Ulb,'ub',Uub);                                            
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% Disturbance  Bounds 
    wub = 0.1;                                    % Upper bound of additive noise value
    wlb = -0.1;                                   % Lower bound of additive noise value

end
