
load BKKmirrorss.mat;
load BKKoriginals.mat;
devFs = 44100;

prompt0 = 'Participant Nr: ';
subid = input(prompt0);

    
    
for t=1:20

trial_wavfile_name = sprintf('BKKmirrorss_pt%d_trial%d.wav', subid, t);

audiowrite(trial_wavfile_name, transpose(BKKmirrors{1,t}), devFs, 'BitsPerSample', 16);

%audiowrite(tap_wavfile_name, transpose(tapdata), devFs, 'BitsPerSample', 16);

end