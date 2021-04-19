%% Evaluate a required binomial sum for bounds
% Monimoy Bujarbaruah

function [val] = eval_binBound(AnomD, N_thres, N, delmatD, Fxrow, boldAvbar, nx)

    val = 0; 
    for j = 1:N-1
        for i = max(N_thres, j):N-1
            val = val + nchoosek(i,j)*norm(AnomD,inf)^(i-j)*norm(delmatD, inf)^j*...
                                    norm(Fxrow*boldAvbar(:,(i-1)*nx*N+1: i*nx*N),1);
        end
    end

end
