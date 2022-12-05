function [sol] = CG(A,b)
    %Solve for Ax = b using Conjugate gradient method
    if nargin~=2
        error ("2 inputs required");
    elseif length(b)~=size(A,1)
        error ("Size mismatch");
    elseif A~=A'
        error ("Not symmetric");
    else
        if iscolumn(b)==0
            b = b';
        end
        tol = 1e-6;     %Tollerence level
        sol = b;        %Initial guess for soln
        res = b - A*sol;%Residue to minimise
        p = res;        %Initial search direction
        %c = 0;
        while (1)
            v = A*p;    %Matrix multiplication
            a = (res'*res)/(p'*v);  %Dot product
            sol = sol + a*p;    %Update soln
            new_res = res - a*v;    %Update residue
            g = (new_res'*new_res)/(res'*res);
            p = new_res + g*p;  %Update search direction
            res = new_res;
            %c = c+1;
            if norm(res)<tol
                break;
            end
        end
        %disp(c);
    end
end

