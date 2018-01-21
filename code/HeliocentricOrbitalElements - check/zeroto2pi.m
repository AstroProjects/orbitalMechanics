function y = zeroto2pi( x )
if x < 0
    y=x+2*pi;
else
    while x>2*pi
        x=x-2*pi;
    end
    y=x;
end
end

