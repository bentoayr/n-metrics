function K = knlPQ2K(KP, KQ, Eg1, Eg2, n1, n2, m1, m2, nn)

I11 = repmat(Eg1(1, :)', 1, m2);
I12 = repmat(Eg1(2, :)', 1, m2);
I21 = repmat(Eg2(1, :), m1, 1);
I22 = repmat(Eg2(2, :), m1, 1);
I1 = sub2ind([n1 n2], I11(:), I21(:));
I2 = sub2ind([n1 n2], I12(:), I22(:));
idx1 = I1(:);
idx2 = I2(:);
vals = KQ(:);

idx1 = [idx1; (1 : nn)'];
idx2 = [idx2; (1 : nn)'];
vals = [vals; KP(:)];

K = sparse(idx1, idx2, vals, nn, nn);
