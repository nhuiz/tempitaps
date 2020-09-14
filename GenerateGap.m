function soundout = GenerateGap(sf, soundin, gapdur, location)
%
% soundout = GenerateEnvelope(sf, soundin, gapdur)
%
% This function applies a single to-be-detected gap with raised cosine 
% gates to a sound. By default, gates have identical duration of 5-ms. 
%
% SF: sample frequency in Hz of the input sound
% SOUNDIN: an array containing the input sound
% GAPDUR: the duration of the gap in ms. 
% LOCATION: the time of gap onset in ms relative to the beginning of the
%           sound

% % cut one sample in case the length of the sound array is an odd number
% if rem(length(soundin), 2)~=0 && size(soundin, 2)==1
%     soundin = soundin(1:end-1);
% elseif rem(length(soundin), 2)~=0 && size(soundin, 2)==2
%     soundin = soundin(1:end-1, 1:2);
% end

rampdur = 3;

rampdurationinsamples = round(rampdur*(sf/1000));
gapdurationinsamples = round(gapdur*(sf/1000));
locationinsamples = round(location*(sf/1000));
% if (onsetdurationinsamples+offsetdurationinsamples) > length(soundin)
%     error ('The duration of onset and offset gates exceedes the overall duration of the sound');
% end

gapramp = (sin(pi*(3/2):pi/(rampdurationinsamples-1):pi*(5/2))+1)/2;
gap = zeros(1,gapdurationinsamples);
sustain = zeros(1, gapdurationinsamples-(2 * length(gapramp)));
if (gapdurationinsamples == 0)
    gap = [fliplr(gapramp),gapramp];
    if (length(gap) > gapdurationinsamples)
        rampdurationinsamples = gapdurationinsamples / 2;
        gapramp = (sin(pi*(3/2):pi/(rampdurationinsamples-1):pi*(5/2))+1)/2;
        gap = zeros(1,2*length(rampdurationinsamples));
        gap = [fliplr(gapramp),gapramp];
    end
else
    gap = [fliplr(gapramp), sustain, gapramp];
end
    
envelope = [ones(1,locationinsamples),gap,ones(1,(length(soundin)-(locationinsamples + length(gap))))];
%plot(envelope)
% apply the gates
if size(soundin, 2)==2
    soundout = [envelope' .* soundin(:, 1), envelope' .* soundin(:, 2)];
else
    if (size(soundin,1) == 1)
        soundout = envelope .* soundin;
    else
        soundout = envelope'.* soundin;
    end
end