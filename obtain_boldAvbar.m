%% Forming boldAvbar Matrix
    
    function  boldAvbar = obtain_boldAvbar(N, nx)
    
    % boldAvbar = [A_1, A_2, ..., A_{N-1}] in (N*nx X N(N-1)nx) , 
    if N >=2
    
        mat = zeros(nx*N, nx*N, N-1); 
       
        for n = 1:N-1         % need from A_1 to A_{N-1}            
            for j =1:N   
                if ((j-1)*nx + n*nx +1)<=nx*N
                    mat((j-1)*nx + n*nx +1: j*nx + n*nx, (j-1)*nx +1: j*nx, n) = eye(nx);     
                end
            end          
        end
    
        boldAvbar = [];
        for k = 1:N-1
            boldAvbar = [boldAvbar, mat(:,:,k)];
        end
    
    else
        boldAvbar = zeros(nx,nx);       % For N=1 case 
         
   
    end