function K = conKnlGphKU(KP, KQ, Eg1, Eg2)

[n1, n2] = size(KP);
[m1, m2] = size(KQ);
nn = n1 * n2;

K = knlPQ2K(KP, KQ, Eg1, Eg2, n1, n2, m1, m2, nn);