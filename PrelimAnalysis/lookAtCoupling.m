addpath ~/gitOther/tensor_toolbox/
bases = RCT.modelBuilder.setupBasis();
B = bases.spkHist.B;

T = size(B,2);

for ff = 1:3
    fname = sprintf("PrelimAnalysis/s%d.mat", ff);
    if(exist(fname, "file"))
        load(fname, "samples");
        S = numel(samples.sample_idx);
        P = size(samples.W,1);
        
        U1 = size(samples.Groups(4).T{2},1);
        U2 = size(samples.Groups(5).T{2},1);
        
        
        coupling_r = zeros(T,  P, P, S);
        coupling_1 = zeros(T, U1, P, S);
        coupling_2 = zeros(T, U2, P, S);
        
        spkHist = zeros(T, P, S);
        
        for ss = 1:S
            if(ss == 1 || ss == S || mod(ss,10) == 0)
                fprintf("sample %d / %d\n", ss, S);
            end
            coupling_r(:,:,:,ss) = double(ktensor(ones(32,1), {samples.Groups(3).T{1}(:,:,ss), samples.Groups(3).T{2}(:,:,ss), samples.Groups(3).V(:,:,ss)}));
            coupling_1(:,:,:,ss) = double(ktensor(ones(16,1), {samples.Groups(4).T{1}(:,:,ss), samples.Groups(4).T{2}(:,:,ss), samples.Groups(4).V(:,:,ss)}));
            coupling_2(:,:,:,ss) = double(ktensor(ones(16,1), {samples.Groups(5).T{1}(:,:,ss), samples.Groups(5).T{2}(:,:,ss), samples.Groups(5).V(:,:,ss)}));
        
            spkHist(:,:,ss) = samples.B(:,:,ss);
        end
        for pp = 1:P
            spkHist(:,pp,:) = spkHist(:,pp,:) + reshape(coupling_r(:,pp,pp,:), [T, 1, S]) ;
            coupling_r(:,pp,pp,:) = 0; 
        end
        fprintf("saving\n");
        save(sprintf("PrelimAnalysis/f%d.mat", ff), "-v7.3", "spkHist", "coupling_r", "coupling_1", "coupling_2", "B");
        fprintf("done\n");
    end
end