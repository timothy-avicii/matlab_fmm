
clear
load test_wo_ans.mat

plot(kernel(2001,:),"b","LineWidth",2)
energy_range_list = [] ;
for i = 1 : 4001
    if kernel(2001, i) < 1e-9
        for j = 1 : 5
            energy_range_list(j, i) = i;
        end
    elseif kernel(2001, i) < 1e-8
        for j = 1 : 4
            energy_range_list(j, i) = i;
        end
    elseif kernel(2001, i) < 1e-7
        for j = 1 : 3
            energy_range_list(j, i) = i;
        end
    elseif kernel(2001, i) < 1e-6
        for j = 1 : 2
            energy_range_list(j, i) = i;
        end
    elseif kernel(2001 ,i) < 1e-5
        energy_range_list(1, i) = i;
    end
end
