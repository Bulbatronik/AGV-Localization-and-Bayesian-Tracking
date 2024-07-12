function h = measurementModel( UE , AP)
numberOfAP = 6;
h = zeros(1, numberOfAP-1);

%d2 = sqrt( sum([UE-AP(2,:)].^2 , 2 ) );

t = 2;

for a = 1:numberOfAP
    if a ==1
        %h(1, a) = sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) - sqrt( sum([UE-AP(2,:)].^2 , 2 ) );
        h(1, a) = -(sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) - sqrt( sum([UE-AP(2,:)].^2 , 2 ) ));
    elseif a > 2
        %h(1, t) = sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) - sqrt( sum([UE-AP(2,:)].^2 , 2 ) );
        h(1, t) = -(sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) - sqrt( sum([UE-AP(2,:)].^2 , 2 ) ));
        t=t+1;
    end
    %i = id(a);
           
end

