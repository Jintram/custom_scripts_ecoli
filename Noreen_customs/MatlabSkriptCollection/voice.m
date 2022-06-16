function beer = ishetalvrijdag

NET.addAssembly('System.Speech');
S = System.Speech.Synthesis.SpeechSynthesizer;
S.Rate = -1;

[d,day] = weekday(now,'long');
D = clock;

hour = D(4);
if hour > 12
    hour = hour - 12;
    ampm = 'PM';
else
    ampm = 'AM';
end

minute = D(5);


str = sprintf('It is %s, %d past %d %s and time to think about beer',day,minute,hour,ampm);

S.Speak(str)

end