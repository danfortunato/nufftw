ntrials = 1000;
ns = 2.^(2:21);
ts = [];

for n = ns
    tic;
    for i = 1:ntrials
        fftw('dwisdom',[]);
        fftw('swisdom',[]);
        in = complex(zeros(n,1));
        out = fft(in);
    end
    tend = toc;
    ts = [ts tend];
    fprintf('# n=%d: %gs\n', n, tend/ntrials);
end
