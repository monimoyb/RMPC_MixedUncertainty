%% Evaluate the chunks of cost for the binomial series sum in the boundings

function [cost02, cost03m] = eval_costChunks(Anom, N_thres, N, Fxrow, boldAvbar, nx, setdelA, delAv)

        cost02 = 0;
        for i = N_thres:N-1
           cost02 = cost02 + norm(kron(eye(N),Anom^i),inf)*norm(Fxrow*boldAvbar(:,(i-1)*nx*N+1: i*nx*N),1);  
        end        
        
        cost03 = zeros(size(delAv,2)/nx, 1); 
        %%% chunks based on deltaA now 
        for i = 1:size(delAv,2)/nx
            cost03(i) = eval_binBound(kron(eye(N),Anom), N_thres, N, kron(eye(N),setdelA(:,:,i)), Fxrow, boldAvbar, nx); 
        end
        cost03m = max(cost03); 
        
end