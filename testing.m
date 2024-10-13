load test_wo_ans.mat

%count = nnz(mask)
count = 0
    for i = 1: size(mask,1)
        for j = 1: size(mask,2)
            if mask(i,j) == 1
                count = count+1;
            end
        end
    end

count