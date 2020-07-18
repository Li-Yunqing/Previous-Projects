#include <stdio.h>
#include <math.h>

#define MAX 100
#define X 10
const double accu_single = 0.0000001;
const double accu_double = 0.000000000000001;


void Bessel_gen_no_rescaling_fp32(float x, float* Jn_fp32, int M);
void Bessel_gen_with_rescaling_fp32(float x, float* Jn_fp32, int M);
void Bessel_gen_no_rescaling_fp64(double x, double* Jn_fp64, int M);
void Bessel_gen_with_rescaling_fp64(double x,double* Jn_fp64, int M);



int main(){
    int i, c;
    float Jn_single_no_rscl[MAX];
    float Jn_single_with_rscl[MAX];
    double Jn_double_no_rscl[MAX];
    double Jn_double_with_rscl[MAX];
    float x_single = X;
    double x_double = X;
    int M = 18;

    Bessel_gen_no_rescaling_fp32(x_single, Jn_single_no_rscl, M);
    Bessel_gen_with_rescaling_fp32(x_single, Jn_single_with_rscl, M);
    Bessel_gen_no_rescaling_fp64(x_double, Jn_double_no_rscl, M);
    Bessel_gen_with_rescaling_fp64(x_double, Jn_double_with_rscl, M);

    //for( c = 0; c < M; c++ ) printf( "REFERENCE!   First  %d      _jn( %d, x )  %.16f\n", c, c, _jn( c, X ) );
    for (i=0;i<M;i++) printf("%dth term of Jn_fp32 without rescaling: %.16lf  \n", i, Jn_single_no_rscl[i]);
    for (i=0;i<M;i++) printf("%dth term of Jn_fp32 with rescaling: %.16lf  \n", i, Jn_single_with_rscl[i]);
    for (i=0;i<M;i++) printf("%dth term of Jn_fp64 without rescaling: %.16lf  \n", i, Jn_double_no_rscl[i]);
    for (i=0;i<M;i++) printf("%dth term of Jn_fp64 with rescaling: %.16lf  \n", i, Jn_double_with_rscl[i]);

    getchar();
    return 0;
}



void Bessel_gen_no_rescaling_fp32(float x, float* Jn_fp32, int M){
    int i;
    float norm_coeff;
    float temp_sum = 0;

    //Backwards calculating of Jn
    Jn_fp32[M] = 0;
    Jn_fp32[M-1] = 1;
    for (i=M; i>=2; i--){
        Jn_fp32[i-2] = 2*(i-1)/x * Jn_fp32[i-1] - Jn_fp32[i];
    }
    
    //Normalize
    for (i=0; i<M; i++){
        temp_sum = temp_sum + 2* Jn_fp32[i]*Jn_fp32[i];
    }
    temp_sum = temp_sum - Jn_fp32[0] * Jn_fp32[0];
    norm_coeff = 1 / sqrt(temp_sum);
    for (i=0; i<M; i++) Jn_fp32[i] = Jn_fp32[i] * norm_coeff;

}


void Bessel_gen_with_rescaling_fp32(float x, float* Jn_fp32, int M){
    int i,k;
    float norm_coeff;
    float temp_sum = 0;
    float rescale_coeff;

    //Backwards calculating of Jn
    Jn_fp32[M] = 0;
    Jn_fp32[M-1] = 1;
    for (i=M; i>=2; i--){
        Jn_fp32[i-2] = 2*(i-1)/x * Jn_fp32[i-1] - Jn_fp32[i];

        //if the new term gets out of dynamic range, rescale.
        if(Jn_fp32[i-2]* accu_single > 1){
            rescale_coeff = 1/(Jn_fp32[i-2]* accu_single);
            for (k=i-2; k<=M; k++) Jn_fp32[k] = Jn_fp32[k] * rescale_coeff;
        }
    }
    //Normalize
    for (i=0; i<M; i++){
        temp_sum = temp_sum + 2* Jn_fp32[i]*Jn_fp32[i];
    }
    temp_sum = temp_sum - Jn_fp32[0] * Jn_fp32[0];
    norm_coeff = 1 / sqrt(temp_sum);
    for (i=0; i<M; i++) Jn_fp32[i] = Jn_fp32[i] * norm_coeff;

}

void Bessel_gen_with_rescaling_fp64(double x, double* Jn_fp64, int M){
    int i,k;
    double norm_coeff;
    double rescale_coeff;
    double temp_sum = 0;

    //Backwards calculating of Jn
    Jn_fp64[M] = 0;
    Jn_fp64[M-1] = 1;
    for (i=M; i>=2; i--){
        Jn_fp64[i-2] = 2*(i-1)/x * Jn_fp64[i-1] - Jn_fp64[i];
        
        //if the new term gets out of dynamic range, rescale.
        if(Jn_fp64[i-2]* accu_double > 1){
            rescale_coeff = 1/(Jn_fp64[i-2]* accu_double);
            for (k=i-2; k<=M; k++) Jn_fp64[k] = Jn_fp64[k] * rescale_coeff;
        }
    }
    //Normalize
    for (i=0; i<M; i++){
        temp_sum += 2* Jn_fp64[i]*Jn_fp64[i];
    }
    temp_sum -= Jn_fp64[0] * Jn_fp64[0];
    norm_coeff = 1.0 / sqrt(temp_sum);
    for (i=0; i<M; i++) Jn_fp64[i] = Jn_fp64[i] * norm_coeff;

}

void Bessel_gen_no_rescaling_fp64(double x, double* Jn_fp64, int M){
    int i,k;
    double norm_coeff;
    double rescale_coeff;
    double temp_sum = 0;

    //Backwards calculating of Jn
    Jn_fp64[M] = 0;
    Jn_fp64[M-1] = 1;
    for (i=M; i>=2; i--){
        Jn_fp64[i-2] = 2*(i-1)/x * Jn_fp64[i-1] - Jn_fp64[i];
        }

    //Normalize
    for (i=0; i<M; i++){
        temp_sum += 2* Jn_fp64[i]*Jn_fp64[i];
    }
    temp_sum -= Jn_fp64[0] * Jn_fp64[0];
    norm_coeff = 1.0 / sqrt(temp_sum);
    for (i=0; i<M; i++) Jn_fp64[i] = Jn_fp64[i] * norm_coeff;

}

