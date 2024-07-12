function h = measurementModel( parameters , UE , AP , TYPE )
h = zeros(1, parameters.numberOfAP);

if strcmp(TYPE,'TDOA')
    h = zeros(1, parameters.numberOfAP-1);
end

for a = 1:parameters.numberOfAP
   switch TYPE
       case 'TOA'
           h(1, a) = sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) ;
       case 'AOA'
           h(1, a) = atan2( ( UE(2)-AP(a,2) ) , ( UE(1)-AP(a,1) ) );
       case 'RSS'
            h(1, a) = parameters.Pref - 10*parameters.np*log10( sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) );
        case 'TDOA'
            if a>1
                h(1, a-1) = sqrt( sum([ UE - AP(a,:) ].^2 , 2 ) ) - sqrt(sum([UE-AP(1,:)].^2 , 2 ) );
            end
   end
end

