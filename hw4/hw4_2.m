threeInd = zeros();
ind4 = 1;
fourInd = zeros();
ind5 = 1;
fiveInd = zeros();
ind6 = 1;
sixInd = zeros();
ind7 = 1;
sevenInd = zeros();
ind8 = 1;
eightInd = zeros();
ind9 = 1;
nineInd = zeros();

for k=1:n
    if(unique_groups(group_idx(k)) == 0)
        zeroInd(ind0) = k;
        ind0 = ind0 + 1;
    elseif(unique_groups(group_idx(k)) == 1)
        oneInd(ind1) = k;
        ind1 = ind1 + 1;
    elseif(unique_groups(group_idx(k)) == 2)
        twoInd(ind2) = k;
        ind2 = ind2 + 1;
    elseif(unique_groups(group_idx(k)) == 3)
        threeInd(ind3) = k;
        ind3 = ind3 + 1;
    elseif(unique_groups(group_idx(k)) == 4)
        fourInd(ind4) = k;
        ind4 = ind4 + 1;
    elseif(unique_groups(group_idx(k)) == 5)
        fiveInd(ind5) = k;
        ind5 = ind5 + 1;
    elseif(unique_groups(group_idx(k)) == 6)
        sixInd(ind6) = k;
        ind6 = ind6 + 1;
    elseif(unique_groups(group_idx(k)) == 7)
        sevenInd(ind7) = k;
        ind7 = ind7 + 1;
    elseif(unique_groups(group_idx(k)) == 8)
        eightInd(ind8) = k;
        ind8 = ind8 + 1;
    elseif(unique_groups(group_idx(k)) == 9)
        nineInd(ind9) = k;
        ind9 = ind9 + 1;
    end
end

%%
zs = fullproj(1:iters,zeroInd);
ones = fullproj(1:iters,oneInd);
twos = fullproj(1:iters,twoInd);
threes = fullproj(1:iters,threeInd);
fours = fullproj(1:iters,fourInd);
fives = fullproj(1:iters,fiveInd);
sixs = fullproj(1:iters,sixInd);
sevens = fullproj(1:iters,sevenInd);
eights = fullproj(1:iters,eightInd);
nines = fullproj(1:iters,nineInd);


%%
[images, labels] = mnist_parse('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');



%%
[images1, labels1] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');
