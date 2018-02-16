for y =1:9
    for x=1:8
        %i = (y-1)*8 + x
        rgb(y, x, 3) = FGS1(y,x);
        rgb(y, x, 2) = FMS(y,x);
        rgb(y, x, 1) = FGS2(y,x);
        
    end
end

image(rgb)