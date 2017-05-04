function getMach(areaRatio)
	% Converging section 
	mach_tolaerence = 1;
	MachConvergingLower = 0;
	MachConvergingUpper = 1;
	MachConverging = 0.5*(MachConvergingLower+MachConvergingUpper);
	G = 1.4;
	% MachConverging This will always be below 1 so initially let's keep it 100
	% calculating the diverging Mach first
	while(mach_tolaerence > 0.001)
		MachConverging = 0.5*(MachConvergingUpper + MachConvergingLower);
		area_difference = areaRatio - ((G+1)/2)^(-((G+1)/(2*(G-1)))) *...
		 ((1+0.5*(G-1)*MachConverging*MachConverging))^(((G+1)/(2*(G-1))))/MachConverging;
		if(area_difference>0)
			MachConvergingUpper = MachConverging;
		else if(area_difference<=0)
			MachConvergingLower = MachConverging;
		mach_tolaerence = MachConvergingUpper - MachConvergingLower;		
		end
	end
end
MachConverging