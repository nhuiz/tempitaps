function soundout2 = linramp3(x,dur,sf)

dur=.075;
f=1000;

t = 1/sf:1/sf:dur;
ramp = linspace(0,1,length(t));

soundout2 = x .* [ramp, ones(1,length(x) - length(ramp))];
soundout2 = soundout2 .* [ones(1,length(x) - length(ramp)), fliplr(ramp)];

t=1/sf:1/sf:dur;

sin(2.*3.14.*f.*t)
end