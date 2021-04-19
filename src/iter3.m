function [x, j, g, s] = iter3(M, b, x0, acc, omega)
% iter3 - iterative methods for tridiagonal matrices
    [~, j] = jacobi(M, b, x0, acc);
    [~, g] = gaussSeidel(M, b, x0, acc);
    [x, s] = sor(M, b, x0, acc, omega);
end


function [x, numIters] = jacobi(M, b, x0, acc)
    x = x0;
    numIters = 0;
    
    while maxNorm(M, x, b) >= acc
        xShiftLeft = circshift(x, 1);
        xShiftRight = circshift(x, -1);
        
        x = (b - M(:, 1) .* xShiftLeft - M(:, 3) .* xShiftRight) ./ M(:, 2);

        numIters = numIters + 1;
    end
end


function [x, numIters] = gaussSeidel(M, b, x0, acc)
    x = x0;
    numIters = 0;
    [len, ~] = size(x);
    
    while maxNorm(M, x, b) >= acc
        x(1) = (b(1) - M(1, 3) * x(2)) ./ M(1, 2);
        for i = 2 : len - 1
            x(i) = (b(i) - M(i, 1) * x(i - 1) - M(i, 3) * x(i + 1)) ./ M(i, 2);
        end
        x(end) = (b(end) - M(end, 1) * x(end - 1)) ./ M(end, 2);
        
        numIters = numIters + 1;
    end
end

function [x, numIters] = sor(M, b, x0, acc, omega)
    x = x0;
    numIters = 0;
    [len, ~] = size(x);
    
    while maxNorm(M, x, b) >= acc
        x(1) = (omega / M(1, 2)) * (b(1) - M(1, 3) * x(2)) + (1 - omega) * x(1);
        for i = 2 : len - 1
            x(i) = (omega / M(i, 2)) * (b(i) - M(i, 1) * x(i - 1) - M(i, 3) * x(i + 1)) + (1 - omega) * x(i);
        end
        x(end) = (omega / M(end, 2)) * (b(end) - M(end, 1) * x(end - 1)) + (1 - omega) * x(1);
        
        numIters = numIters + 1;
    end
end


function r = maxNorm(M, x, b)
    xShiftLeft = circshift(x, 1);
    xShiftRight = circshift(x, -1);
    
    calcB = M(:, 1) .* xShiftLeft + M(:, 2) .* x + M(:, 3) .* xShiftRight;
    
    r = max(abs(calcB - b));
end