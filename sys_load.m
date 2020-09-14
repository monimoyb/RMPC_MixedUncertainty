%% Function defining System and Constraint Matrices 


function [Anom,Bnom, delAv, delBv, K, A, B, X, U, Xlb, Xub, Ulb, Uub, nx, nu, wub,wlb, Q, R, N_max, N_thres] = sys_load(vertflag, epsA, epsB)

%% Considering two states and one scalar input 

Anom = [1, 0.15; 0.1, 1];
Bnom = [0.1; 1.1];                    
nx = size(Anom,2); nu = size(Bnom,2); 

%%% Now read the vertflag. Value 1 means full set of vertices necessary, given a norm
%%% representation of uncertainty. Value 0 means structured uncertainty representing vertices.  

if vertflag == 1
    %%%%%%%%% Get the error vertices matrices
     delAv = [ [epsA, 0; epsA, 0], [epsA, 0; 0 epsA], [epsA, 0; -epsA, 0], [epsA, 0; 0, -epsA],...
                 [0, epsA; epsA, 0], [0, epsA; 0, epsA], [0, epsA; -epsA, 0], [0, epsA; 0, -epsA],...
                 [-epsA, 0; epsA, 0], [-epsA, 0; 0 epsA], [-epsA, 0; -epsA, 0], [-epsA, 0; 0, -epsA],...
                 [0, -epsA; epsA, 0], [0, -epsA; 0, epsA], [0, -epsA; -epsA, 0], [0, -epsA; 0, -epsA]];    
                    
    delBv = [[0; -epsB], [0; epsB], [epsB; 0], [-epsB; 0] ];
     
     %%%%%%%%% choose the feedback gain K here%%%
     K = place(Anom, Bnom, [0.72; 0.65]);           
     
else
     delAv = [ [0, epsA; epsA, 0], [0, epsA; -epsA, 0], [0, -epsA; epsA, 0], [0, -epsA; -epsA, 0]];     
            
     delBv = [[0; -epsB], [0; epsB], [epsB; 0], [-epsB; 0] ];
       
     %%%%%%%%% choose the feedback gain K here%%%
     K = place(Anom, Bnom, [0.735; 0.75]);             
end
                                                                           

%% Set the true A and B matrices (satisfy the above bounds)
A = [1, 0.1; 
       0, 1.0 ]; 

B = [0; 
     1.0];                                                                                                                       
%% Weights
Q =  10*eye(nx);
R =   2*eye(nu);
%% Horizon and cut-off for bounds
N_max = 4;            % try horizon from 5 to 1 by shrinking
N_thres = 4;           % beyond this value do the binomial tail. 
 
%% Considering constraints of the form -a<=x(i)<=a and -ulb<=u<=uub
Xlb = -[8; 8];
Xub = -Xlb; 
Ulb =  -4; 
Uub = -Ulb; 

X = Polyhedron('lb',Xlb,'ub',Xub);
U = Polyhedron('lb',Ulb,'ub',Uub);                                            
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Disturbance  Bounds 
wub = 0.1;                                   % Upper bound of additive noise value
wlb = -0.1;                                   % Lower bound of additive noise value

end