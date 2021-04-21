%% Compute all the bounding terms. First go up to threshold horizon and then bound binomial tail 
% Monimoy Bujarbaruah

function [t_w, t_1, t_2, t_3, t_delTaA, t_delTaB] = bounds(Fx, Anom, Bnom, N, N_thres, boldAvbar, delAv, delBv, nx, nu)

    t_w = zeros(size(Fx,1),1);
    t_0 = zeros(size(Fx,1),1);
    t_1 = zeros(size(Fx,1),1); 
    t_2 = zeros(size(Fx,1),1); 
    t_3 = zeros(size(Fx,1),1); 
    t_delTaA = zeros(size(Fx,1),1); 
    t_delTaB = zeros(size(Fx,1),1); 

    %% extract the set of possible error matrices vertex wise %%%% 
    setdelA = zeros(nx, nx, size(delAv,2)/nx);
    setdelB = zeros(nx, nu, size(delBv,2)/nu);

    for i = 1:size(delAv,2)/nx
        setdelA(:,:,i) = delAv(:,(i-1)*nx + 1: i*nx);  
    end

    for i = 1:size(delBv,2)/nu
        setdelB(:,:,i) = delBv(:,(i-1)*nu + 1: i*nu);  
    end

    %% Now extract the set of all possible matrices A_\Delta for bounding
    setA = zeros(nx, nx, size(delAv,2)/nx);
    for i = 1:size(delAv,2)/nx
        setA(:,:,i) = Anom + delAv(:,(i-1)*nx + 1: i*nx);  
    end

    %% These are needed later    
    costdelB = zeros(size(delBv,2)/nu,1); 
    for j = 1: size(delBv,2)/nu
        costdelB(j) = norm(kron(eye(N), setdelB(:,:,j)),inf);
    end
    t_delb = max(costdelB); 

    costdelA = zeros(size(delAv,2)/nx,1); 
    for j = 1: size(delAv,2)/nx
        costdelA(j) = norm(setdelA(:,:,j),inf);
    end
    t_dela = max(costdelA); 


    boldA1bar = zeros(nx*N, nx*N);
    for j = 1:N
        tmpMat = [];
        for k = 1:j
            tmpMat = [tmpMat, Anom^(j-k)];
        end
        boldA1bar(nx*(j-1)+1: nx*j, 1:nx*j) = tmpMat;    
    end

    for row = 1: size(Fx,1)
        for i = 1:size(delAv,2)/nx
            tmp(i) = norm(Fx(row,:)*boldA1bar*kron(eye(N), setdelA(:,:,i)), 1); 
        end
        t_delTaA(row) = max(tmp); 

        for i = 1:size(delBv,2)/nu
           tmp(i) = t_delb*norm(Fx(row,:)*boldA1bar*kron(eye(N), setdelB(:,:,i)), 1); 
        end   
        t_delTaB(row) = max(tmp); 
    end


    %% Start all cases here on

    if N>=2 && N<=N_thres         % give values of N \geq 2 and N \leq N_thres. 

        %% FORM POSSIBLE VERTICES FOR A_DELA^n from n=1 to N-1
        ApowVertMatrCell{1} = setA; 
        for i = 2: N-1
            matind  = permn([1:size(delAv,2)/nx], i);
            ApowVertMatr = zeros(nx,nx,size(matind,1));
            for j = 1: size(matind,1)
                ApowVertMatr(:,:,j) = eye(nx); 
                for k = 1:size(matind,2)
                    ApowVertMatr(:,:,j) = ApowVertMatr(:,:,j)*setA(:,:,matind(j,k));
                end
            end
            ApowVertMatrCell{i} = ApowVertMatr;
        end

        %%%%%%%%%% ALL VERTEX MATRICES STORED %%%%%
        %% Now store the stacked matrices for all combinations  
        % find all possible combinations of vertex indices for each power to
        % N-1. Max will occur at one of these combinations   
        for i = 1:N-1
           tmp = ApowVertMatrCell{i};
           sz{i} = [1:size(tmp,3)];
        end
        %%% hard coding the line below!!!%%%%%
        % need to do it from N-1 onward to 1.         
        if N == 4
            pre_final = combvec(sz{3}, sz{2}, sz{1})';                                                       % N = 4 case
        elseif N == 3
            pre_final=combvec(sz{2}, sz{1})';                                                                % N = 3 case
        elseif N == 2
            pre_final=combvec(sz{1})';                                                                       % N = 2 case 
        end

        for i=1:N-1
            final_combvert(:,i)=pre_final(:,(N-1)-i+1);
        end

        %%% have all the combinations now. Need to go to each cell array and pick out the corresponding matrices
        % by looping over rows of the final_combvert matrix!!! 
        mat_comb = zeros(nx, (N-1)*nx, size(final_combvert,1) ); 
        for j = 1:size(final_combvert,1) 
           tmp_comb = [];
           for k = 1: (N-1)
                tmp_comb = [ tmp_comb, ApowVertMatrCell{k}(:,:,final_combvert(j,k))]; 
           end
           mat_comb(:,:,j) = tmp_comb; 
        end
        % Forming all the stacked matrices here for all the combinations
        delmat = zeros(N*nx*(N-1), N*nx,size(final_combvert,1)); 
        for j = 1: size(final_combvert,1)
            delmattmp = []; 
            for k = 1:N-1
                delmattmp = [delmattmp; kron( eye(N), mat_comb(:,(k-1)*nx+1: k*nx,j) - Anom^k )];               
            end
            delmat(:,:,j) = delmattmp; 
        end
        %%% this one is for tw
        delmat_tw = zeros(N*nx*(N-1), N*nx,size(final_combvert,1)); 
        for j = 1: size(final_combvert,1)
            delmattmp = []; 
            for k = 1:N-1
                delmattmp = [delmattmp; kron( eye(N), mat_comb(:,(k-1)*nx+1: k*nx,j) )];
            end
            delmat_tw(:,:,j) = delmattmp; 
        end

        %% Find the bounds by trying ALL vertex combinations row-wise
        for row = 1: size(Fx,1)                                                                                    
            % formed stacked matrices. Each one is a combination 
            % Now need to unpack and evaluate cost !!!
            cost_w = zeros(size(final_combvert,1), 1); 
            for j = 1: size(final_combvert,1)
                 cost_w(j) = norm(Fx(row,:)*boldAvbar*delmat_tw(:,:,j), 1);
            end    
            %%%%%%%%%%% finding t_w %%%%%%%%%%%    
            t_w(row) = max(cost_w); 

            % form the rest from t_0 to t_3 now %%%
            cost0 = zeros(size(final_combvert,1), 1); 
            for j = 1: size(final_combvert,1)
                 cost0(j) = norm(Fx(row,:)*boldAvbar*delmat(:,:,j), 1);
            end    
            %%%%%%%%%%% finding t_0 %%%%%%%%%%%    
            t_0(row) = max(cost0); 
            %%% Now finding t_1, t_2 %%%
            t_1(row) = t_0(row)*t_dela;
            t_2(row) = t_0(row)*t_delb;       
            %%% this is used for t_3  
            cost3 = zeros(size(final_combvert,1), 1); 
            for j = 1: size(final_combvert,1)
                 cost3(j) = norm(Fx(row,:)*boldAvbar*delmat(:,:,j)*kron(eye(N), Bnom), 1);
            end
            %%% finding t_3
            t_3(row) = max(cost3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

    elseif N==1    % For N=1 case   

          disp('**N=1 case encountered. Nothing to bound with norms in this case**')

    else           % N>N_thres. Use binomila tail for the remaining bounds
        %% FORM POSSIBLE VERTICES FOR A_DELA^n from n=1 to N_thres-1 
        ApowVertMatrCell{1} = setA; 
        for i = 2: N_thres-1
            matind  = permn([1:size(delAv,2)/nx], i);
            ApowVertMatr = zeros(nx,nx,size(matind,1));
            for j = 1: size(matind,1)
                ApowVertMatr(:,:,j) = eye(nx); 
                for k = 1:size(matind,2)
                    ApowVertMatr(:,:,j) = ApowVertMatr(:,:,j)*setA(:,:,matind(j,k));
                end
            end
            ApowVertMatrCell{i} = ApowVertMatr;
        end

        %%%%%%%%%% ALL VERTEX MATRICES STORED %%%%%
        %% Now store the stacked matrices for all combinations  
        % find all possible combinations of vertex indices for each power to
        % N_thres-1. Max will occur at one of these combinations   
        for i = 1:N_thres-1
           tmp = ApowVertMatrCell{i};
           sz{i} = [1:size(tmp,3)];
        end
        %%% hard coding the line below!!!%%%%%
        % need to do it from N_thres-1 onward to 1. 
        if N_thres == 4
            pre_final=combvec(sz{3}, sz{2}, sz{1})';                                                       % N_thres = 4 case
        elseif N_thres == 3
            pre_final=combvec(sz{2}, sz{1})';                                                              % N_thres = 3 case
        elseif N_thres == 2
            pre_final=combvec(sz{1})';                                                                     % N_thres = 2 case 
        end

        for i=1:N_thres-1
            final_combvert(:,i)=pre_final(:,(N_thres-1)-i+1);
        end

        %%% have all the combinations now. Need to go to each cell array and pick out the corresponding matrices
        % by looping over rows of the final_combvert matrix!!! 
        mat_comb = zeros(nx, (N_thres-1)*nx, size(final_combvert,1) ); 
        for j = 1:size(final_combvert,1) 
           tmp_comb = [];
           for k = 1: (N_thres-1)
                tmp_comb = [ tmp_comb, ApowVertMatrCell{k}(:,:,final_combvert(j,k))]; 
           end
           mat_comb(:,:,j) = tmp_comb; 
        end

        % Forming all the stacked matrices here for all the combinations
        delmat = zeros(N*nx*(N_thres-1), N*nx,size(final_combvert,1)); 
        for j = 1: size(final_combvert,1)
            delmattmp = []; 
            for k = 1:N_thres-1
                delmattmp = [delmattmp; kron( eye(N), mat_comb(:,(k-1)*nx+1: k*nx,j) - Anom^k)];             
            end
            delmat(:,:,j) = delmattmp; 
        end

        % This one is for t_w
        delmat_tw = zeros(N*nx*(N_thres-1), N*nx,size(final_combvert,1)); 
        for j = 1: size(final_combvert,1)
            delmattmp = []; 
            for k = 1:N_thres-1
                delmattmp = [delmattmp; kron( eye(N), mat_comb(:,(k-1)*nx+1: k*nx,j) )];
            end
            delmat_tw(:,:,j) = delmattmp; 
        end

        %% Find the bounds by trying ALL vertex combinations row-wise
        for row = 1: size(Fx,1)                                
            % formed stacked matrices. Each one is a combination 
            % Now need to unpack and evaluate cost !!!
            % Form the required chunk of the boldAvbar matrix

            boldAvbar_thres = boldAvbar(:, 1:(N_thres-1)*nx*N); 
            %%% first chunk of cost_w
            costw1 = zeros(size(final_combvert,1), 1); 
            for j = 1: size(final_combvert,1)
                 costw1(j) = norm(Fx(row,:)*boldAvbar_thres*delmat_tw(:,:,j), 1);
            end  

            %%% second chunk of costw from Biomial sum series
            [costAnomGeom, costMixed] = eval_costChunks(Anom, N_thres, N, Fx(row,:), boldAvbar, nx, setdelA, delAv);

            %%%% finding t_w %%%%%%%%%%%%%%%%%
            t_w(row) = max(costw1) + costAnomGeom + costMixed;

            % form the rest needed for t_0 to t_3

            %%% first chunk of cost0
            cost01 = zeros(size(final_combvert,1), 1); 
            for j = 1: size(final_combvert,1)
                 cost01(j) = norm(Fx(row,:)*boldAvbar_thres*delmat(:,:,j), 1);
            end  

            %%% second chunk of cost0 from Biomial sum series %%%%%%%%%%%%%
            Fxrow = Fx(row, :); 
            sum = 0; 
            for j = N_thres+1:N                     % first Nthres sums amount to zero 
                bin_term = 0; 
                for i = 1:j-N_thres
                    for k = 1:j-i
                        bin_term = bin_term + nchoosek(j-i,k)*norm(Anom,inf)^(j-i-k)*t_dela^k;
                    end
                end
                sum = sum + norm(Fxrow(1, (j-1)*nx+1: j*nx), 1)*bin_term;             
            end

            %%%% finding t_0 %%%%%%%%%%%%%%%%%
            t_0(row) = max(cost01) + sum; 

            %%%% finding t_1, t_2, t_3 %%%%%%%%%%%%%%%%%
            t_1(row) = t_0(row)*t_dela;  
            t_2(row) = t_0(row)*t_delb;                                         
            t_3(row) = t_0(row)*norm(kron(eye(N), Bnom), inf);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

    end

end
