a = rand(3,1000);
c = zeros(2,10);
for i = 1:10;
    m = 10*i^2;
    a = rand(3,m);
    Q = rand(3);
    tic
    b = sum(a'*Q.*a',2);
    c(1,i) = toc;
    tic
    b1 = 0.5*a'*Q*a;
    b2 = b1(logical(eye(m)))';
    c(2,i) = toc;
    clear b b1 b2
end
