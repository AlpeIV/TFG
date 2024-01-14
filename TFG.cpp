#include "stdio.h"
#include <cmath>
#include <stdlib.h>

//Constants

double Ergs_In_Joule = 10000000.000000;
double ElectronVolts_in_Erg = 624150907446.000000;
double Pi_Const = 3.141592653589793238462643383279502884197169399375105820974944;
double SecondsInAyear = 31557600.0;
double BoltzmannConstant = 1.380649 * 1e-16; // cm2 g s-2 K-1
double BoltzmannConstant_SI = 1.380649 * 1e-23; // 
double BoltzmannConstant_eV = 8.6173 * 1e-5; // eV K-1
double ElectronMass = 9.1095 * 1e-28; // g
double LightSpeed = 2.9979 * 1e+10; // cm/s
double ElectronMassEnergy = 510998.9460999999 * (1/ElectronVolts_in_Erg); //erg
double ElectronMassEnergy_eV = 510998.9460999999; //eV
double MinimumEtaEffValue = 1.0000000000;
double HydrogenMass = 1.673534 * 1e-24; //g




//Relevant Values

double DistanceBetween_AcceleratorStar = 3.0 * 1e+12;   // cm
double StarBrightness = 1e+39;   // erg/s
double StarSurfaceTemperature = 40000; //Kelvin
double StarMassLoss = ((1.988478e+33 /*Solar Mass in g*/)/SecondsInAyear) * 1e-6;  // 10^-6 Solar Masses/year in g/s
double WindSpeed = 200000000; // cm/s
double AcceleratorEnergyInjection = 1e+36;   // erg/s
double MaxElectronEnergy = 160.2176565000016251 ; // erg
double MinElectronEnergy = 000.001602176565000016251 ; // erg



//Expressions:

double surfaceArea(double Radius){    //cm^2
	return 4.0 * Pi_Const * Radius * Radius;
}

double BrightnessDensity(double Surface, double Brightness){
	return Brightness/Surface;
}

double WindDensity(double urface, double Windspeed, double massloss){
	return massloss/(urface * Windspeed);
}

double CharacteristicNrgy(double Temperatoora){
	return 3.0 * BoltzmannConstant * Temperatoora;
}

double ParticleDeensity(double luminosity, double surfeis, double charact_nrg){
	return luminosity/(surfeis * charact_nrg);
}

bool IsKleinNishina_Regime(double StarSurfaceTemperature_Kelvin, double ElectronEnergy_ergs){
	
	double GammaRelativisticParameter = (ElectronEnergy_ergs/ElectronMassEnergy);
	//printf("Gamma = %.2e \n", GammaRelativisticParameter);
	double KleinNishinaLowerThr_eV = 1.0/((StarSurfaceTemperature_Kelvin * BoltzmannConstant_eV)/ (ElectronMassEnergy_eV));
	//printf("1/kT = %.2e \n", KleinNishinaLowerThr_eV);
	if(KleinNishinaLowerThr_eV < (4.49 * GammaRelativisticParameter)){
		//printf("KN Regime\n");
		return true;
	}else{
		//printf("Th Regime\n");
		return false;
	}
	

	return false;
}


double accelerator_BField_Gauss(double U_star, double Lambdad){

	return sqrt(U_star * Lambdad * 8.0 * Pi_Const);
}

double acceleratorCaractTyme(double EtaEfficiency, double ParticleEnergy_TeV, double BFieldAccelerator_Gauss){
	if(EtaEfficiency < 1.00000000000000){
		printf("Eta cannot be less than 1, fixing...\n");
		EtaEfficiency = 1.00000000000000;
	}
	
	//printf("Particle Energy:\n%f TeV\n", ParticleEnergy_TeV);
	//printf("Magnetic Field:\n%.2e G\n", BFieldAccelerator_Gauss);
	return (0.1 * EtaEfficiency * ParticleEnergy_TeV)/BFieldAccelerator_Gauss;
}



double SynchroCoolingTym(double ParticleEnergy_TeV, double BFieldAccelerator_Gauss){

	return (400.0)/(BFieldAccelerator_Gauss* BFieldAccelerator_Gauss * ParticleEnergy_TeV);
}





double CoolingTym_IC_KN(double U_star, double T_Star, double E_Particles_TeV){
	
	double W_100 = U_star/100.0;
	
	return (100.0/W_100) * pow(T_Star/30000, 1.7000000000000) * pow(E_Particles_TeV, 0.7000000000000);
}



double IC_Th_EdotCalcul(double U_Star, double E_Particles){
	return 0.039* U_Star * E_Particles * E_Particles;  //fabs of -0.039
}

double IC_Th_EdotCalcul_TeVs(double U_Star, double E_Particles_TeV){
	return 0.039* U_Star * E_Particles_TeV * E_Particles_TeV;  //fabs of -0.039
}




double Edot_BremsstrahlungCalc(double n_windParicleAmount, double E_Particle){
	return 1e-15*n_windParicleAmount*E_Particle;
}




double Tyme_fromEdotE(double E_Particles, double Edot){
	return E_Particles/Edot;
}



double findSmallestValue(double a, double b, double c) {
    if (a <= b && a <= c) {
        return a;
    } else if (b <= a && b <= c) {
        return b;
    } else {
        return c;
    }
}

double findSmallestValueDos(double a, double b) {
    if (a <= b) {
        return a;
    } 
    return b;
}




int main(){
	
	long long Counter = 0;
	
	 double SummFracsSy = 0.0;
	double SummFracsIC = 0.0;
	double SummFracsBr = 0.0;
	double IComptonMinimum = 100000000000000.0, IComptonEenergy;
	
	
	bool TriggeredCoolEq = false;
	
	printf("Distance\n%.6e cm\n\n", DistanceBetween_AcceleratorStar);
		
	double Surface = surfaceArea(DistanceBetween_AcceleratorStar);
	printf("Surface Area\n%.2e cm^2\n\n", Surface);
	
	printf("Star Luminosity\n%.2e erg / s\n\n", StarBrightness);
	
	double LuminosityAtAccelerator  = BrightnessDensity(Surface, StarBrightness);
	printf("Luminosity Flux at Interest Region\n%.2e erg / s cm^2 \n\n", LuminosityAtAccelerator);
	double U_star = LuminosityAtAccelerator/LightSpeed;
	printf("Energy Density U* at Interest Region\n%.2e erg / cm^3 \n\n", U_star);
	
	
	printf("Star Mass Loss M'\n%.2e g / s \n\n", StarMassLoss);
	printf("Star Wind Speed\n%.2e cm / s \n\n", WindSpeed);
	
	double WindDensityValue  = WindDensity(Surface,  WindSpeed, StarMassLoss );
	printf("Wind Density\n%.2e g / cm^3 \n\n", WindDensityValue);
		
	double CharacteristicEnergyStar  = CharacteristicNrgy(StarSurfaceTemperature);
	printf("Star Characteristic Energy\n%.2e erg \n", CharacteristicEnergyStar);
	printf("%f eV \n\n", CharacteristicEnergyStar * ElectronVolts_in_Erg);
	
	double EscapeTime = DistanceBetween_AcceleratorStar/LightSpeed;
	printf("Escape Time %.2e s \n\n",  EscapeTime);
	
	
	printf("\n\n\nWorking with Electron Particles\n\n");

	
	double LambdaEquivalenceRateMinimum = 0.01;
	double LambdaEquivalenceRateMax = 100;
	double LambdaEquivalenceRate = LambdaEquivalenceRateMinimum ;
	double WorkingEta = MinimumEtaEffValue ;
	printf("Assuming Accelerator Efficiency of:\n%f \n\n", WorkingEta);
	
	//printf("Max Electron Energy\n%.2e erg \n", ElectronEnergy);
	//printf("%f TeV \n\n", ElectronEnergy * ElectronVolts_in_Erg * 1e-12);
	
	double ElectronEnergy = MinElectronEnergy ;
	

	
	
	FILE* AllResults = fopen("DataTable.dat", "w");
	if(AllResults == NULL){
		printf("Error in File Write");
		return 1;
	}
	

	
	//fprintf(AllResults, "Electron Energy     Characteristic Time      IC Time        Synchro Time        Bremss       Escape\n");	
	
	
	do{
		
			printf("\n\nU_B/U_*\n%.2f\n", LambdaEquivalenceRate);
			fprintf(AllResults, "Lambda relation U_B/U_*     %06.2f\n\n", LambdaEquivalenceRate);
			fprintf(AllResults, "Electron Energy     Characteristic Accelerator Time      IC Time        Synchro Time        Bremss       Escape\n");
			double ICTime_ForPrint = 0;
	
	
		do{
		
			//printf("Electron Energy\n%.2e erg \n", ElectronEnergy);
			double ElectronEnergy_TeV =  ElectronEnergy * ElectronVolts_in_Erg * 1e-12;
			//printf("%f TeV\n", ElectronEnergy_TeV);
	
			
			double MagneticFieldAtAccelerator = accelerator_BField_Gauss(U_star, LambdaEquivalenceRate);
			//printf("Magnetic Field\n%.2e Gauss \n", MagneticFieldAtAccelerator);
	
			double AcceleratorCharacteristicTime = acceleratorCaractTyme(WorkingEta, ElectronEnergy_TeV, MagneticFieldAtAccelerator );
			//printf("Characteristic Accelerator Time \n%.2e s \n", AcceleratorCharacteristicTime);
		
		
		
			//Synchrotron
			double SynchroCoolTime = SynchroCoolingTym(ElectronEnergy_TeV, MagneticFieldAtAccelerator);
			//printf("Synchrotron Cooling Time \n%.2e s \n", SynchroCoolTime);
			if(SynchroCoolTime > AcceleratorCharacteristicTime){
				//printf("Synchrotron Energy Loss is faster than gain, secondary process\n");
			}else{
				//printf("Synchrotron will be a primary process\n");
			}
			//double SynchroMaxParticleNrgy = SynchroMaxEnergy( WorkingEta,  MagneticFieldAtAccelerator );
			//printf("Synchrotron Maximum Particle Energy \n%f TeV \n", SynchroMaxParticleNrgy);
			
			
		    bool WeInKN = IsKleinNishina_Regime(StarSurfaceTemperature, /*ergs*/ ElectronEnergy);
			if(WeInKN){
				//printf("Inverse Compton in KN\n");
			
				double ICTime = CoolingTym_IC_KN(U_star, StarSurfaceTemperature, ElectronEnergy_TeV);
				ICTime_ForPrint = ICTime;
				//printf("IC Cooling Time KN \n%.2e s \n", ICTime);
				//fprintf(IC_TimeResults, "%.2e\n", ICTime);
			
				if(ICTime > AcceleratorCharacteristicTime){
					//printf("IC Energy Loss is faster than gain, secondary process\n");
				}else{
					//printf("IC Emission will be predominant\n");
				}
				

		
			}else{
				//printf("Inverse Compton in Th\n");
				double Edot = IC_Th_EdotCalcul( U_star,  ElectronEnergy);
				double IC_Th_CoolinTime = Tyme_fromEdotE(ElectronEnergy, Edot);
				ICTime_ForPrint = IC_Th_CoolinTime;
				//printf("IC Cooling Time Th \n%.2e s \n", IC_Th_CoolinTime);
				//fprintf(IC_TimeResults, "%.2e\n", IC_Th_CoolinTime);
			}
			
			
			//Bremsstrahlung 
			double Edot_Bremsstrahlung = Edot_BremsstrahlungCalc( WindDensityValue / HydrogenMass, /*ergs*/ ElectronEnergy);
			//printf("Bremsstrahlung Edot \n%.2e erg/s \n", Edot_Bremsstrahlung);
			double TimeCharcteristic_Bremsstrahlung = Tyme_fromEdotE(ElectronEnergy, Edot_Bremsstrahlung);
			//printf("Bremsstrahlung Time \n%.2e s \n", TimeCharcteristic_Bremsstrahlung);
		
	
			//printf("\n");
			
			//File Print
			fprintf(AllResults, "%.2e          %.2e                           %.2e      %.2e           %.2e    %.2e\n", ElectronEnergy, (AcceleratorCharacteristicTime), (ICTime_ForPrint), (SynchroCoolTime), (TimeCharcteristic_Bremsstrahlung), EscapeTime);
			
			//Logic for fraction computation
			
			
			//Computing luminosity fraction correspondence
			
			double TimeEfficient = findSmallestValue(SynchroCoolTime, ICTime_ForPrint, TimeCharcteristic_Bremsstrahlung); 
			Counter++;
			
			
			//Synchrotron
			double SigmaSy;
			if(SynchroCoolTime < EscapeTime){
				SigmaSy = 1.0;
			}else{
				SigmaSy = EscapeTime/SynchroCoolTime;
			}
			SummFracsSy = SummFracsSy + (SigmaSy*(TimeEfficient/SynchroCoolTime));

			
			//Inverse Compton
			double SigmaIC;
			if(ICTime_ForPrint < EscapeTime){
				SigmaIC = 1.0;
			}else{
				SigmaIC = EscapeTime/ICTime_ForPrint;
			}
			SummFracsIC = SummFracsIC + (SigmaIC*(TimeEfficient/ICTime_ForPrint)); 
	
	
			//Bremsstrahlung
			double SigmaBr;
			if(TimeCharcteristic_Bremsstrahlung < EscapeTime){
				SigmaBr = 1.0;
			}else{
				SigmaBr = EscapeTime/TimeCharcteristic_Bremsstrahlung;
			}
			SummFracsBr = SummFracsBr + (SigmaBr*(TimeEfficient/TimeCharcteristic_Bremsstrahlung)); 
			
			
			
            //Logic for energy calculations in t_eq
            double TimeCool = findSmallestValue(SynchroCoolTime, ICTime_ForPrint, TimeCharcteristic_Bremsstrahlung);
            if((TimeCool <= AcceleratorCharacteristicTime) && (TriggeredCoolEq == false)){
                TriggeredCoolEq = true;
                printf("Equality found for Lambda = %6.2f,  at Electron Energy %.2e,   Cool Time %.2e\n", LambdaEquivalenceRate, ElectronEnergy, TimeCool);
    
                //printf("B: %.2e\n", MagneticFieldAtAccelerator);
            }
			
			if(IComptonMinimum > ICTime_ForPrint){
				IComptonMinimum = ICTime_ForPrint;
				IComptonEenergy = ElectronEnergy;
			}

			
			ElectronEnergy = ElectronEnergy + MinElectronEnergy;
			if(ElectronEnergy > (MaxElectronEnergy + 0.001)){		//Maximum Limit
				break;
			}
		}while(ElectronEnergy < (MaxElectronEnergy + 0.001));   //To test other possible small regions
		
		
		printf("IC Most efficient time %.2e s  at %.2e erg\n", IComptonMinimum, IComptonEenergy);
		
		SummFracsSy = 100 * SummFracsSy/Counter;
		printf("Luminosity fraction for Sy %6.4f\%\n", SummFracsSy);
		
		SummFracsIC = 100 * SummFracsIC/Counter;
		printf("Luminosity fraction for IC %6.4f\%\n", SummFracsIC);
		
		SummFracsBr = 100 * SummFracsBr/Counter;
		printf("Luminosity fraction for Br %.2e\%\n", SummFracsBr);
		

		TriggeredCoolEq = false;
		printf("\n\n");
		
		Counter = 0;
		
		SummFracsSy = 0.0;
		SummFracsIC = 0.0;
		SummFracsBr = 0.0;
		IComptonMinimum = 100000000000000000000.0;;
		
		
		//Reset to beginning of energy range
		ElectronEnergy = MinElectronEnergy;
		
		//Increment Lambda
		LambdaEquivalenceRate = LambdaEquivalenceRate * 10;
		if(LambdaEquivalenceRate > LambdaEquivalenceRateMax){  //Maximum Limit
			break;
		}
			
	}while(LambdaEquivalenceRate < (LambdaEquivalenceRateMax + 3));  //To test other possible smaller values
	

	fclose(AllResults);

	return 0;
}
