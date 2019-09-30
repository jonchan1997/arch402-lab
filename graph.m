n = [64,128, 256, 512, 1024];
title('N vs. Gflops')
xlabel('N') 
ylabel('Gflops')
legend({'UNOPTIMIZED_DGEMM','SSE_DGEMM'},'Location','southwest')
unop_gflops = [0.679970, 0.506629, 0.513142, 0.726387, 0.304151];
sse_gflops = [2.319645, 2.028151, 2.390608, 3.041345, 1.106510];
plot(n,unop_gflops,n,sse_gflops);