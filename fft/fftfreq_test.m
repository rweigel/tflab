assert(isempty(fftfreq(0)))
assert(all(fftfreq(1) == 0))
assert(all(fftfreq(2) == [0,-1/2]))
assert(all(fftfreq(3) == [0,1/3,-1/3]))
assert(all(fftfreq(4) == [0,1/4,-1/2,-1/4]))
assert(all(fftfreq(5) == [0,1/5,2/5,-2/5,-1/5]))

f = fftfreq(10);
assert(all(f(f>=0) == fftfreqp(10)))