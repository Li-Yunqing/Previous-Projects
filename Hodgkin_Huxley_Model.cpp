/*

	Details to specify: m n h V初始值;
			    输入;

*/

#include <iostream>
#include <cmath>

#define E 2.71828182845904523536
	#define MembraneCapitance 1
	#define G_Na 120
	#define G_K 36
	#define G_L 0.33
	#define V_Na 50
	#define V_K -77
	#define V_L -54.4



using namespace std;


void TEMPARATURE_DEPENDENT_FUNCTION(const double temparature,double * temparatureFactor){
	*temparatureFactor = pow(3,((temparature-6.3)/10));
}


void calculateMDerivative(double m,double V,double * md,double temparatureFactor) {
	*md = temparatureFactor * 0.1 * ((25-V) / (pow(E,((25-V)/10)-1))) *(1-m)
		- temparatureFactor * 4 * pow(E,(-V/18));

}

void calculateNDerivative(double n,double V,double * nd, double temparatureFactor) {
	*nd = temparatureFactor * 0.01 * ((10-V) / (pow(E,((10-V)/10)-1))) *(1-n)
		- temparatureFactor * 0.125 * pow(E,(-V/80));
}

void calculateHDerivative(double h,double V,double * hd, double temparatureFactor) {
	*hd = temparatureFactor * 0.07 * (pow(E,(-V/20))) *(1-h)
		- temparatureFactor * (1/((pow(E,((30-V)/10))))+1);
}

void calculateVDerivative(double inputFunction, double V,double m,double n,double h,double * Vd,double temparatureFactor) {
	*Vd = - G_Na * pow(m,3) * h * (V-(V_Na))
		- G_K * pow(n,4) * (V-V_K)
		- G_L * (V-V_L)
		+ inputFunction;
}

void calculateNewValue(double * currentValue,double Derivative,double stepLength) {
	*currentValue += Derivative * stepLength;
}



int main() {

	double V,Vd,m,md,n,nd,h,hd;
	double* referenceToV=&V;
	double* referenceToVd=&Vd;
	double* referenceToM=&m;
	double* referenceToMd=&md;
	double* referenceToN=&n;
	double* referenceToNd=&nd;
	double* referenceToH=&h;
	double* referenceToHd=&hd;
	
	V=-40; m=1; n=1; h=1;

	double timeLength,numberOfSteps;
	double temparature,temparatureFactor;
	double* refToTemparatureFac = &temparatureFactor;

	
	double inputFunction = 3;

	cin >> timeLength >> numberOfSteps;
	cin >> temparature;
	const double stepLength = timeLength/numberOfSteps;

	TEMPARATURE_DEPENDENT_FUNCTION(temparature,refToTemparatureFac);

	

	for (int i=0; i!=numberOfSteps; i++) {
		calculateMDerivative(m,V,referenceToMd,temparatureFactor);
		calculateNDerivative(n,V,referenceToNd,temparatureFactor);
		calculateHDerivative(h,V,referenceToHd,temparatureFactor);
		calculateVDerivative(inputFunction,V,m,n,h,referenceToVd,temparatureFactor);

		calculateNewValue(referenceToM,md,stepLength);
		calculateNewValue(referenceToN,nd,stepLength);
		calculateNewValue(referenceToH,hd,stepLength);
		calculateNewValue(referenceToV,Vd,stepLength);

		cout << "The current value of V is " << V << endl;

}
	

return 0;
}

