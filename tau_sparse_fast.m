function [Tsum,ketfinal] = tau_sparse_fast(r,Ti,structU,ketold,tMax)

    found = false;
    tit  = length(Ti);
    Tsum = 0;
    Uold = speye(length(ketold));
    it = 0;
    tol = 1e-3;
    ketfinal = ketold;

    while not(found)

        ketnew = structU(tit,1).matrix*ketold;

        %mod = ketnew'*ketnew;
        mod = norm(ketnew);
        if abs(mod-r)<tol
            ketfinal = ketnew;
            Tsum = Tsum + Ti(tit);
            break;
        end

        if Tsum > tMax
            break;
        end
        
        if mod < r
            tit = tit-1;
            if tit == 0
                found = true;
            end
        else
            ketold = ketnew;
            Tsum = Tsum + Ti(tit);
        end
        it = it + 1;
    end
    
    ketfinal = Uold*ketold;
    
end

