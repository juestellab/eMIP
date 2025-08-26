function M = extendND(M, N)
    if nargin < 2
        N = 3;
    end
    [sz1, sz2] = size(M);
    if sz2 == 1 && sz1 > 1
        M = M';
    end
    M(:,sz2+1:N) = M(:,sz2);
end

