function  [P, score, un_rounded_score] = g_align_distance(all_graphs, n, k, nuc_vs_psd, should_round)

%nuc norm 1

    A = zeros(n,n,k);
    for i = 1:k
        A(:,:,i) = all_graphs.adj{i};
    end

%     A = zeros(n,n,k);
%     for i = 1:k
%         A(:,:,i) = all_graphs{i};
%     end

    cvx_begin quiet

    variable P(n*k,n*k)

    s = 0;
    for i = 1:k
        for j = 1:k
            s = s + norm(A(:,:,i)*P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) ) - P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) )*A(:,:,j) , 'fro');
        end
    end

    minimize (0.5*s   )

    subject to  
            if (nuc_vs_psd == 1)
                norm_nuc(P) <= n*k;
            else
            	P == semidefinite(n*k);    
            end
            
            diag(P) == 1;
            for i = 1:k
                for j = 1:k
                    P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) ) >= 0;
                    sum(P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) )) == 1;
                    sum(P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) )') == 1;
                end
            end
    cvx_end

    un_rounded_score = cvx_optval/(k*k);
    
    full_FP = P;
    
    if (should_round == 1)

        for i = 1:k-1
            for j = i+1:k
                cvx_begin quiet
                    variable FP(n,n);
                    minimize (  -trace(FP*((P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) ) )' + 0.001*rand(n,n))) )

                    subject to
                        FP*ones(n,1) == 1;
                        FP'*ones(n,1) == 1;
                        0 <= FP;
                        FP <= 1;
                cvx_end
                full_FP([1:n] +  n*(i-1) , [1:n] +  n*(j-1) ) = round(FP);
                full_FP([1:n] +  n*(j-1) , [1:n] +  n*(i-1) ) = round(FP)';
            end
        end

    end
    
    P = full_FP;
    
    s = 0;
    for i = 1:k
        for j = 1:k
            s = s + norm(A(:,:,i)*P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) ) - P([1:n] +  n*(i-1) , [1:n] +  n*(j-1) )*A(:,:,j) , 'fro');
        end
    end
    
    score = 0.5*s;
    
end