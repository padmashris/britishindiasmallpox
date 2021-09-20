load solarMFmagnitudes
fb = cwtfilterbank('SignalLength',length(sm),'SamplingPeriod',hours(1));
scaleSpectrum(fb,sm)