/*
Author: Albert Patterson
email: apatterson@alumni.cmu.edu


This code features a demonstration of a class used to simulate two-port RF
circuits. With this class, a user can build circuits from passive components
and simulate the S parameters in order to characterize the frequency response
of the circuit. Furthermore, the center frequency and insertion loss of a
filter design may be calculated.

Included here is the class as well as a demonstration of a simulation of a
third-order filter composed of two-port contour-mode resonators. The user will
prompted to give a quality factor value, which will be used in the simulation
of the filter. A typical quality factor is ~1500. The user will be given the
option of saving the simulated S parameters of |S21| values in tab-delineated
files, which can easily be opened in excel to visualize the frequency response
of the circuit.
*/


#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

/*
This class is used for modeling the frequency response a two port network
consisting of cascaded elements or capacitors, inductors and resistors.

Two port networks can be cascaded together to form complex networks.

Scattering parameter and transmission parameter values can be calculate

Simulation can me run for any frequency range

The insertion loss and center frequency of a filter design can be calculated
*/
class TwoPortNetwork{
    private:
        // type of parameters held, T for transmission and S for scattering
        char param_type_;
        // number of frequency values
        int n_freqs_;
        // frequencies used in simulation in Hz
        double* freqs_;
        // termination impedance is 50 this is a critical parameter in design
        // of RF circuits
        int r0_=50;
        //values of two port network parameters
        double *(param_values_[8]);

        // calculate admittance of capacitor, inductor, or resistor
        double** calcAdmittance(char component_type,float component_val);
        // calculate impedance of capacitor, inductor, or resistor
        double** calcImpedance(char component_type,float component_val);
        // calculate the real part of the product of two complex numbers
        double calcComplexProductReal(double real_val_1,double imag_val_1,
                                      double real_val_2, double imag_val_2);
        // calculate the imaginary part of the product of two complex numbers
        double calcComplexProductImag(double real_val_1,double imag_val_1,
                                      double real_val_2, double imag_val_2);
        // calculate the real part of the quotient of two complex numbers
        double calcComplexQuotientReal(double real_val_1,double imag_val_1,
                                       double real_val_2, double imag_val_2);
        // calculate the imaginary part of the product of two complex numbers
        double calcComplexQuotientImag(double real_val_1,double imag_val_1,
                                       double real_val_2, double imag_val_2);
        // find the maximum value in an array
        double findMax(double* arr,int arr_len);
        // fint the index of the maximum value in an array
        int findMaxIndex(double* arr,int arr_len);
    public:
        //construct and initialize parameters as transmission parameters
        // based on provided circuit information
        TwoPortNetwork(int n_freqs_in,double *freqs_in, char component_type,
                       float component_val,bool is_series);
        //construct by copying another two port network
        TwoPortNetwork(TwoPortNetwork& tpn_a);
        // destroy two port network
        ~TwoPortNetwork();

        //get the type of parameters stored
        const char get_param_type();
        // get a parameter value;
        double get_param_value(int row, int col);
        //convert parameters from transmission to scattering
        void convert_params_to_s();
        //get forward transmission in dB
        double* get_s21_dB();
        // get the insertion loss (il) of the two port network
        double get_insertion_loss();
        // get the center frequency (f0) of the two port network
        double get_center_frequency();
        // connect two port network tpn_b to the existing two port network
        void cascade(TwoPortNetwork& tpn_b);
};



// calculate the admittance value for a resistor inductor or capacitor
double** TwoPortNetwork::calcAdmittance(char component_type,
                                        float component_val){

    // admittance will be a 2 x n_freqs_ array
    // where the first column is the real part
    // and the second column is the imaginary part
    double **admittance=new double*[2];
    admittance[0]=new double[n_freqs_];
    admittance[1]=new double[n_freqs_];

    // check if component is a capacitor
    if (component_type=='c'){
        for(int i=0;i<n_freqs_;i++){
            admittance[0][i]=0;//real
            admittance[1][i]=2*3.14*freqs_[i]*component_val;//imaginary
        }
    // check if component is an inductor
    }else if(component_type=='l'){
        for(int i=0;i<n_freqs_;i++){
            admittance[0][i]=0;//real
            admittance[1][i]=-1/(2*3.14*freqs_[i]*component_val);//imaginary
        }
    // check if component is a resistor
    }else if(component_type=='r'){
        for(int i=0;i<n_freqs_;i++){
            admittance[0][i]=1/component_val;//real
            admittance[1][i]=0;//imaginary
        }
    }else{
        // warn user of error
        cout<<"Component must be char c, l or r for capacitor, inductor"
            <<"or resistor, respectively. Admittance not calculated.";
    }

    return admittance;
}

// calculate the impedance value for a resistor inductor or capacitor
double** TwoPortNetwork::calcImpedance(char component_type,
                                       float component_val){

    // impedance will be a 2 x n_freqs_ array
    // where the first column is the real part
    // and the second column is the imaginary part
    double **impedance=new double*[2];
    impedance[0]=new double[n_freqs_];
    impedance[1]=new double[n_freqs_];

    // check if component is a capacitor
    if (component_type=='c'){
        for(int i=0;i<n_freqs_;i++){
            impedance[0][i]=0;//real
            impedance[1][i]=-1/(2*3.14*freqs_[i]*component_val);//imaginary
        }
    // check if component is an inductor
    }else if(component_type=='l'){
        for(int i=0;i<n_freqs_;i++){
            impedance[0][i]=0;//real
            impedance[1][i]=2*3.14*freqs_[i]*component_val;//imaginary
        }
    // check if component is a resistor
    }else if(component_type=='r'){
        for(int i=0;i<n_freqs_;i++){
            impedance[0][i]=component_val;//real
            impedance[1][i]=0;//imaginary
        }
    }else{
        // warn user of error
        cout<<"Component must be char c, l or r for capacitor, inductor"
            <<"or resistor, respectively. Impedance not calculated.";
    }

    return impedance;
}

// calculate the real part of (r1+j*i1)*(r2+j*i2)=
// (r1*r2-i1*i2)+j*(r1*i2+r2+i1)
// where j=sqrt(-1)
double TwoPortNetwork::calcComplexProductReal(double real_val_1,
                                              double imag_val_1,
                                              double real_val_2,
                                              double imag_val_2){

    double real=real_val_1*real_val_2-imag_val_1*imag_val_2;
    return real;
}

// calculate the imaginary part of
// (r1+j*i1)*(r2+j*i2)=(r1*r2-i1*i2)+j*(r1*i2+r2+i1)
// where j=sqrt(-1)
double TwoPortNetwork::calcComplexProductImag(double real_val_1,
                                              double imag_val_1,
                                              double real_val_2,
                                              double imag_val_2){

    double imag=real_val_1*imag_val_2+imag_val_1*real_val_2;
    return imag;
}

// calculate the real part of
// (r1+j*i1)/(r2+j*i2)=
// ((r1*r2+i1*i2)/(r2^2+i2^2))+j*((r2*i1-r1+i2)/(r2^2+i2^2))
// where j=sqrt(-1) and r1, r2, i1, and i2 are real numbers
double TwoPortNetwork::calcComplexQuotientReal(double real_val_1,
                                               double imag_val_1,
                                               double real_val_2,
                                               double imag_val_2){

    double real=(real_val_1*real_val_2+imag_val_1*imag_val_2)/
                (pow(real_val_2,2)+pow(imag_val_2,2));

    return real;
}

// calculate the imaginary part of
// (r1+j*i1)/(r2+j*i2)=
// ((r1*r2+i1*i2)/(r2^2+i2^2))+j*((r2*i1-r1+i2)/(r2^2+i2^2))
// where j=sqrt(-1) and r1, r2, i1, and i2 are real numbers
double TwoPortNetwork::calcComplexQuotientImag(double real_val_1,
                                               double imag_val_1,
                                               double real_val_2,
                                               double imag_val_2){

    double imag=(real_val_2*imag_val_1-real_val_1*imag_val_2)/
                (pow(real_val_2,2)+pow(imag_val_2,2));

    return imag;
}

// find max element of array
double TwoPortNetwork::findMax(double* arr,int arr_len){
        double max_val=arr[0];
        for(int i=1;i<arr_len;i++){
            if (arr[i]>max_val){
                max_val=arr[i];
            }
        }
        return max_val;
}

// find max element of array
int TwoPortNetwork::findMaxIndex(double* arr,int arr_len){
        int max_ind=0;
        double max_val=arr[0];
        for(int i=1;i<arr_len;i++){
            if (arr[i]>max_val){
                max_val=arr[i];
                max_ind=i;
            }
        }
        return max_ind ;
}



// TwoPortNetwork constructor
// transmission matrix values assigned depending on the frequency and
// component values
// freqs: frequencies in Hz
// component: 'r' for resistor, 'c' for capacitor, 'l' for inductor
// is_series: true for series connection, false for parallel connection
TwoPortNetwork::TwoPortNetwork(int n_freqs_in,double *freqs_in,
                               char component_type,float component_val,
                               bool is_series){
    // preinitialize values
    param_type_='T';
    n_freqs_=n_freqs_in;
    freqs_=new double[n_freqs_];
    for(int row=0;row<n_freqs_;row++){
        freqs_[row]=freqs_in[row];
    }

    for(int col=0;col<8;col++){
        param_values_[col]=new double[n_freqs_];
    }

    // assign values as transmission paremeters
    if (is_series){
        double **z=calcImpedance(component_type,component_val);
        for(int row=0;row<n_freqs_;row++){
            //a11 real
            param_values_[0][row]=1;
            //a11 imaginary
            param_values_[1][row]=0;
            //a12 real
            param_values_[2][row]=z[0][row];
            //a12 imaginary
            param_values_[3][row]=z[1][row];
            //a21 real
            param_values_[4][row]=0;
            //a21 imaginary
            param_values_[5][row]=0;
            //a22 real
            param_values_[6][row]=1;
            //a22 imaginary
            param_values_[7][row]=0;
        }

    }else{
        double **y=calcAdmittance(component_type,component_val);
        for(int row=0;row<n_freqs_;row++){
            //a11 real
            param_values_[0][row]=1;
            //a11 imaginary
            param_values_[1][row]=0;
            //a12 real
            param_values_[2][row]=0;
            //a12 imaginary
            param_values_[3][row]=0;
            //a21 real
            param_values_[4][row]=y[0][row];;
            //a21 imaginary
            param_values_[5][row]=y[1][row];;
            //a22 real
            param_values_[6][row]=1;
            //a22 imaginary
            param_values_[7][row]=0;
        }
    }
}

//construct by copying another two port network
TwoPortNetwork::TwoPortNetwork(TwoPortNetwork& tpn_a){
    param_type_=tpn_a.param_type_;
    n_freqs_=tpn_a.n_freqs_;
    freqs_=new double[n_freqs_];
    for(int row=0;row<n_freqs_;row++){
            freqs_[row]=tpn_a.freqs_[row];
    }

    for(int col=0;col<8;col++){
        param_values_[col]=new double[n_freqs_];
    }

    for(int row=0;row<n_freqs_;row++){
        for(int col=0;col<8;col++){
            param_values_[col][row]=tpn_a.param_values_[col][row];
        }
    }
}

//destroy two port network
TwoPortNetwork::~TwoPortNetwork(){
    delete[] freqs_;
    for(int col=0;col<8;col++){
        delete[] param_values_[col];
    }
    delete[] param_values_;
}



// return the type of parameters held
const char TwoPortNetwork::get_param_type(){
    return param_type_;
}

// get a parameter value;
double TwoPortNetwork::get_param_value(int row, int col){
    return param_values_[col][row];
};



// calculate s parameters from transmission parameters a11, a12, a21, a22
// s11=(a11+a12/r0-a21*r0-a22)/(a11+a12/r0+a21*r0+a22)
// s12=2*(a11*a22- a12*a21)/(a11+a12/r0+a21*r0+a22)
// s21=2/(a11+a12/r0+a21*r0+a22)
// s22=(-a11+a12/r0-a21*r0+a22)/(a11+a12/r0+a21*r0+a22)
void TwoPortNetwork::convert_params_to_s(){
    if(param_type_=='T'){
        //  the new parameter type will be S for scattering
        this->param_type_='S';

        double den_real,den_imag,a11_real,a11_imag,
            a12_real,a12_imag,a21_real,a21_imag,a22_real,a22_imag,
            ref_prod_real,ref_prod_imag,
            tran_prod_real,tran_prod_imag,
            s11_num_real,s11_num_imag,
            s12_num_real,s12_num_imag,
            s21_num_real,s21_num_imag,
            s22_num_real,s22_num_imag,
            s11_real,s11_imag,
            s12_real,s12_imag,
            s21_real,s21_imag,
            s22_real,s22_imag;

        // calculate s parameters for every frequency value
        for(int row=0; row<n_freqs_;row++){
            // store transmission parameter values
            // in easily recognized variables
            a11_real=param_values_[0][row];
            a11_imag=param_values_[1][row];
            a12_real=param_values_[2][row];
            a12_imag=param_values_[3][row];
            a21_real=param_values_[4][row];
            a21_imag=param_values_[5][row];
            a22_real=param_values_[6][row];
            a22_imag=param_values_[7][row];

            // denominator: a11+a12/r0+a21*r0+a22
            // the denominator is the same in all
            // expressions for converting from
            // transmission to scattering parameters
            den_real=a11_real+a12_real/r0_+a21_real*r0_+a22_real;
            den_imag=a11_imag+a12_imag/r0_+a21_imag*r0_+a22_imag;

            // s11 numerator: a11+a12/r0-a21*r0-a22
            s11_num_real=a11_real+a12_real/r0_-a21_real*r0_-a22_real;
            s11_num_imag=a11_imag+a12_imag/r0_-a21_imag*r0_-a22_imag;
            // s11 = s11_num / den
            s11_real=calcComplexQuotientReal(s11_num_real,s11_num_imag,
                                             den_real,den_imag);
            s11_imag=calcComplexQuotientImag(s11_num_real,s11_num_imag,
                                             den_real,den_imag);

            // s12 numerator: 2*(a11*a22- a12*a21)
            // ref_prod=a11*a22
            ref_prod_real=calcComplexProductReal(a11_real,a11_imag,a22_real,
                                                 a22_imag);
            ref_prod_imag=calcComplexProductImag(a11_real,a11_imag,a22_real,
                                                 a22_imag);
            // tran_prod=a12*a21
            tran_prod_real=calcComplexProductReal(a12_real,a12_imag,a21_real,
                                                  a21_imag);
            tran_prod_imag=calcComplexProductImag(a12_real,a12_imag,a21_real,
                                                  a21_imag);
            // s12_num=2*(a11*a22- a12*a21)=2*(ref_prod-tran_prod)
            s12_num_real=2*(ref_prod_real-tran_prod_real);
            s12_num_imag=2*(ref_prod_imag-tran_prod_imag);
            // s12 = s12_num/den
            s12_real=calcComplexQuotientReal(s12_num_real,s12_num_imag,
                                             den_real,den_imag);
            s12_imag=calcComplexQuotientImag(s12_num_real,s12_num_imag,
                                             den_real,den_imag);

            // s21 numerator: 2
            s21_num_real=2;
            s21_num_imag=0;
            //s21 =s21_num/den
            s21_real=calcComplexQuotientReal(s21_num_real,s21_num_imag,
                                             den_real,den_imag);
            s21_imag=calcComplexQuotientImag(s21_num_real,s21_num_imag,
                                             den_real,den_imag);

            // s22 numerator: -a11+a12/r0-a21*r0+a22
            s22_num_real=-a11_real+a12_real/r0_-a21_real*r0_+a22_real;
            s22_num_imag=-a11_imag+a12_imag/r0_-a21_imag*r0_+a22_imag;
            // s22=s22_num/den
            s22_real=calcComplexQuotientReal(s22_num_real,s22_num_imag,
                                             den_real,den_imag);
            s22_imag=calcComplexQuotientImag(s22_num_real,s22_num_imag,
                                             den_real,den_imag);

            // assign s11
            param_values_[0][row]=s11_real;
            param_values_[1][row]=s11_imag;
            //assign s12
            param_values_[2][row]=s21_real;
            param_values_[3][row]=s21_imag;
            //assign s21
            param_values_[4][row]=s12_real;
            param_values_[5][row]=s12_imag;
            //assign s12
            param_values_[6][row]=s22_real;
            param_values_[7][row]=s22_imag;
        }
    }else if(param_type_=='S'){
        //do nothing, we already have scattering parameters
    }else{
        cout<<param_type_;
        perror("Bad parameter type. Only S and T allowed");
        //error
    }
}

// get the magnitude of the forward transmission S21 in dB
double* TwoPortNetwork::get_s21_dB(){

    double *s21db=new double[n_freqs_];
    double s21_real;//real part of forward transmission
    double s21_imag;// imaginary part of forward transmission
    if(param_type_=='S'){//Get forward transmission using S parameters
        for(int row=0; row<n_freqs_;row++){
            //real part of forward transmission
            s21_real=param_values_[4][row];
            // imaginary part of forward transmission
            s21_imag=param_values_[5][row];
            // magnitude of forward transmission in dB
            s21db[row]=20*log10(sqrt(pow(s21_real,2)+pow(s21_imag,2)));
        }
    }else if(param_type_=='T'){//Get forward transmission using T parameters
        this->convert_params_to_s();//Change to S parameters
        s21db=this->get_s21_dB();// get forward transmission
    }

    return s21db;
}

// calculate the insertion loss of a two port network
double TwoPortNetwork::get_insertion_loss(){
        double* s21dbs=this->get_s21_dB();// get |S21|
        double insertion_loss=findMax(s21dbs,n_freqs_);//find max |S21|
        return -1*insertion_loss;// report insertion loss as positive
}

// calculate the center frequency of a two port network
double TwoPortNetwork::get_center_frequency(){
    // get |S21|
    double* s21dbs=this->get_s21_dB();
    //find max |S21|
    int center_frequency_index=findMaxIndex(s21dbs,n_freqs_);
    //find frequency were |S21| is max
    double center_frequency=freqs_[center_frequency_index];
    return center_frequency;
}

// Connnect a two port network tpn_b onto an existing two port network
// Multiply together T matrices to model cascaded two port networks
void TwoPortNetwork::cascade(TwoPortNetwork& tpn_b){
    // temporarily hold new T values of cascaded tpn
    double temp_param_values[8];

    // We must have T parameters to do the calculation
    // for simplicity, let this TwoPortNetwork have T parameters
    // {a11,a12,a21,a22} and let tpn_b have T parameters {b11 b12 b21 b22}
    // where each parameter has a real and an imaginary part
    if ((param_type_='T')&&(tpn_b.get_param_type()=='T')){
        for(int row=0;row<n_freqs_;row++){
            // real part a11*b11+a12*b21
            temp_param_values[0]=calcComplexProductReal(//real part a11*b11
                                        //real part a11
                                        param_values_[0][row],
                                        //imaginary part a11
                                        param_values_[0+1][row],
                                        //real part b11
                                        tpn_b.param_values_[0][row],
                                        //imaginary part b11
                                        tpn_b.param_values_[0+1][row])
                                +calcComplexProductReal(//real part a12*b21
                                        //real part a12
                                        param_values_[2][row],
                                        //imaginary part a12
                                        param_values_[2+1][row],
                                        //real part b21
                                        tpn_b.param_values_[4][row],
                                        //imaginary part b21
                                        tpn_b.param_values_[4+1][row]);
            // imaginary part a11*b11+a12*b21
                                    //imaginary part a11*b11
            temp_param_values[0+1]=calcComplexProductImag(
                                        //real part a11
                                        param_values_[0][row],
                                        //imaginary part a11
                                        param_values_[0+1][row],
                                        //real part b11
                                        tpn_b.param_values_[0][row],
                                        //imaginary part b11
                                        tpn_b.param_values_[0+1][row])
                                //imaginary part a12*b21
                                +calcComplexProductImag(
                                        //real part a12
                                        param_values_[2][row],
                                        //imaginary part a12
                                        param_values_[2+1][row],
                                        //real part b21
                                        tpn_b.param_values_[4][row],
                                        //imaginary part b21
                                        tpn_b.param_values_[4+1][row]);

            //real part of a11*b12+a12*b22
                                    //real part of a11*b12
            temp_param_values[2]=calcComplexProductReal(
                                        //real part of a11
                                        param_values_[0][row],
                                        //imaginary part of a11
										param_values_[0+1][row],
                                        // real part of b12
										tpn_b.param_values_[2][row],
                                        //imaginary part of b12
										tpn_b.param_values_[2+1][row])
                                    //real part of a12*b22
                                +calcComplexProductReal(
                                        // real part of a12
										param_values_[2][row],
                                        // imaginary part of a12
										param_values_[2+1][row],
                                        // real part of b22
										tpn_b.param_values_[6][row],
                                        //imaginary part of b22
										tpn_b.param_values_[6+1][row]);
             //imaginary part of a11*b12+a12*b22
									// imaginary part of a11*b12
			temp_param_values[2+1]=calcComplexProductImag(
                                        // real part of a11
										param_values_[0][row],
                                        // imaginary part of a11
										param_values_[0+1][row],
                                        // real part of b12
										tpn_b.param_values_[2][row],
                                        // imaginary part of b12
										tpn_b.param_values_[2+1][row])
                                // imaginary part of a12*b22
								+calcComplexProductImag(
                                        // real part of a12
										param_values_[2][row],
                                        // imaginary part of a12
										param_values_[2+1][row],
                                        // real part of b22
										tpn_b.param_values_[6][row],
                                        // imaginary part of b22
										tpn_b.param_values_[6+1][row]);

            //real part of a21*b11+a22*b21
                                    // real part of a21*b11
            temp_param_values[4]=calcComplexProductReal(
                                        // real part of a21
										param_values_[4][row],
                                        // imaginary part of a21
										param_values_[4+1][row],
                                        // real part of b11
										tpn_b.param_values_[0][row],
                                        // imaginary part of b11
										tpn_b.param_values_[0+1][row])
                                    // real part of a22*b21
                                +calcComplexProductReal(
                                        // real part of a22
										param_values_[6][row],
                                        // imaginary part of a22
										param_values_[6+1][row],
                                        // real part of b21
										tpn_b.param_values_[4][row],
                                        // imaginary part of b21
										tpn_b.param_values_[4+1][row]);
            //imaginary part of a21*b11+a22*b21
									// imaginary part of a21*b11
			temp_param_values[4+1]=calcComplexProductImag(
                                        // real part of a21
										param_values_[4][row],
                                        // imaginary part of a21
										param_values_[4+1][row],
                                        // real part of b11
										tpn_b.param_values_[0][row],
                                        // imaginary part of b11
										tpn_b.param_values_[0+1][row])
                                // imaginary part of a22*b21
								+calcComplexProductImag(
                                        // real part of a22
										param_values_[6][row],
                                        // imaginary part of a22
										param_values_[6+1][row],
                                        // real part of b21
										tpn_b.param_values_[4][row],
                                        // imaginary part of b21
										tpn_b.param_values_[4+1][row]);

            //real part of a21*b12+a22*b22
                                    // real part of a21*b12
            temp_param_values[6]=calcComplexProductReal(
                                        // real part of a21
										param_values_[4][row],
                                        // imaginary part of a21
										param_values_[4+1][row],
                                        // real part of b12
										tpn_b.param_values_[2][row],
                                        // imaginary part of b12
										tpn_b.param_values_[2+1][row])
                                    // real part of a22*b22
                                +calcComplexProductReal(
                                        // real part of a22
										param_values_[6][row],
                                        // imaginary part of a22
										param_values_[6+1][row],
                                        // real part of b22
										tpn_b.param_values_[6][row],
                                        // imaginary part of b22
										tpn_b.param_values_[6+1][row]);
            //imaginary part of a21*b12+a22*b22
									// imaginary part of a21*b12
			temp_param_values[6+1]=calcComplexProductImag(
                                        // real part of a21
										param_values_[4][row],
                                        // imaginary part of a21
										param_values_[4+1][row],
                                        // real part of b12
										tpn_b.param_values_[2][row],
                                        // imaginary part of b12
										tpn_b.param_values_[2+1][row])
                                // imaginary part of a22*b22
								+calcComplexProductImag(
                                        // real part of a22
										param_values_[6][row],
                                        // imaginary part of a22
										param_values_[6+1][row],
                                        // real part of b22
										tpn_b.param_values_[6][row],
                                        // imaginary part of b22
										tpn_b.param_values_[6+1][row]);

            for(int col=0;col<8;col++){
                param_values_[col][row]=temp_param_values[col];
            }
        }
    }else{
        // Warn the user that parameters must be T and inform that no change
        // was made
        cout<<"Only Transmission parameters can be multiplied together\n";
        cout<<"Elements not connected, no change made to TwoPortNetwork\n";
    }
}






// calculate linearly spaced array of doubles
// from start to finish, returning array of length n_freqs
double* linearlySpacedArray(double start, double finish, int n_freqs);

// request a quality factor value from user
float requestQ();
// check if user wants to save S parameter data
bool saveSData();
// check if user wants to save |S21| data
bool saveS21Data();
// save S parameter data to a file
void saveSDataToFile(double *freqs, TwoPortNetwork &tpn_a,int n_freqs);
// save |S21| data to a file
void saveS21DataToFile(double *freqs, double *s21_vals,int n_freqs);
// Assess the filter quality for the user
void assessDesign(double insertion_loss, double center_frequency);





int main (){

    // setup simulation
    int n_freqs=1000;// number of frequency values in simulation
    double start_freq=0.95e9;//950 MHz
    double end_freq=1.05e9;//1050 MHz

    // get user input
    float Q=requestQ();//What is the quality factor of resonators
    //do you want to save the S parameter data? (all s parameters)
	bool save_s=saveSData();
    bool save_s21=saveS21Data();// do you want to save the S21 data?

    //run simulation
    // frequencies for simulation
	double* freqs=linearlySpacedArray(start_freq, end_freq, n_freqs);
    // motional resistance of resonators (depends on Q)
	float motional_resistance_val=2000/Q*3.43;
    //motional capacitance of resonators
	float motional_capacitance_val=2.32e-14;
    // motional inductance of resonators
	float motional_inductance_val=1.09e-6;
    // static capacitance of resonators
	float static_capacitance_val=3.18e-12;
//
//                 Circuit models simulated in this example
//
//                 Resonator
// Rm=motional resistor, Cm=motional capacitor, Lm=motional inductor,
//                        C0=static capacitor
//
//                    Rm    Cm    Lm
//        o----------^V^V^--||--uuuu---------o
//               _|_                 _|_
// port 1     C0 ___              C0 ___     port 2
//                |                   |
//        o----------------------------------o
//
//                                              Third-Order Filter
//
//                   Rm    Cm   Lm               Rm    Cm   Lm               Rm    Cm   Lm
//        o---------^V^V^--||--uuuu-------------^V^V^--||--uuuu-------------^V^V^--||--uuuu---------o
//              _|_                 _|_     _|_                 _|_     _|_                 _|_
// port 1    C0 ___              C0 ___  C0 ___              C0 ___  C0 ___              C0 ___     port2
//               |                   |       |                   |       |                   |
//        o------------------------------------------------------------------------------------------o

    // create two port networks for each component
    //two port network for the static capacitance
	TwoPortNetwork static_capacitor(n_freqs, freqs,'c',
                                    static_capacitance_val,false);
    //two port network for the motional resistance
	TwoPortNetwork motional_resistor(n_freqs, freqs,'r',
                                    motional_resistance_val,true);
    //two port network for the motional inductance
	TwoPortNetwork motional_inductor(n_freqs, freqs,'l',
                                    motional_inductance_val,true);
    //two port network for the motional capacitance
	TwoPortNetwork motional_capacitor(n_freqs, freqs,'c',
									motional_capacitance_val,true);

    // build a single resonator
    //start out with a static capacitor
	TwoPortNetwork resonator(static_capacitor);
    // connect the motional resistor
	resonator.cascade(motional_resistor);
    // connect the motional inductor
	resonator.cascade(motional_inductor);
    // connect the motional capacitor
	resonator.cascade(motional_capacitor);
    // connect the static capacitor to finish the resonator
	resonator.cascade(static_capacitor);

    //build the filter
    TwoPortNetwork filter(resonator);// start out with a resonator
    filter.cascade(resonator);// connect another resonator
    filter.cascade(resonator);// connect a third resonator

    //calculate the insertion loss of filter
    double insertion_loss=filter.get_insertion_loss();
    // calculate the center frequency of filter
    double center_frequency=filter.get_center_frequency();

    // tell the user about filter performance and if it is ok
    assessDesign(insertion_loss, center_frequency);

    // save full s parameter data
    if (save_s){
        saveSDataToFile(freqs,filter,n_freqs);
    }

    // save s21 data
    if (save_s21){
        double* s21_vals=filter.get_s21_dB();
        saveS21DataToFile(freqs,s21_vals,n_freqs);
    }

    cout<<"\n\n";
    return 0;
}





// calculate linearly spaced array of doubles
// from start to finish and the number of values
// given by n_freqs
double* linearlySpacedArray(double start, double finish, int n_freqs){
    double *arr=new double[n_freqs];// linearly spaced array
    double delt=(finish-start)/(n_freqs-1);//spacing between values
    for(int i=0;i<n_freqs;i++){
        arr[i]=start+delt*i;
    }
    return arr;
}

// request a quality factor value from user
float requestQ(){
    // Request a quality factor value from the user
    cout<<"What is the quality factor of your contour-mode resontors?\nQ=";
    float Q;//quality factor
    cin>>Q;//get quality factor value
    while (Q<=0){
        cout<<"Quality factor must be positive!\nQ=";
        cin>>Q;
    }
    return Q;
}

// check if user wants to save S parameter data
bool saveSData(){
    // Ask if user wants to say S parameters
    cout<<"\nWould you like to save S parameters in text file";
	cout<<"sdata_albert_patterson_code_sample.txt?\n(Y/N): ";
    string question_1_answer;
    cin>>question_1_answer;// get answer from user

    // boolean value for saving s parameter data default false
    bool save_s=false;

    // action take depending on answer
    // answer no
    if((question_1_answer=="n")||(question_1_answer=="no")||
		(question_1_answer=="N")||(question_1_answer=="No")||
		(question_1_answer=="NO")){

        cout<<"\nNot saving data.\n";//user answers no, dont save

    //answer yes
    }else if((question_1_answer=="y")||(question_1_answer=="yes")||
			(question_1_answer=="Y")||(question_1_answer=="Yes")||
			(question_1_answer=="YES")){

        cout<<"\nSaving data.\n";// user answers yes
        save_s=true;// save data

    }else{
        //unintelligible answer, don't save
		cout<<"\nChoice not understood. Not Saving.\n";
    }
    return save_s;
}

// check if user wants to save |S21| data
bool saveS21Data(){
    // Ask if user wants to say S21 data
    cout<<"\nWould you like to save |S21| values in dB in text file";
	cout<<"s21data_albert_patterson_code_sample.txt?\n(Y/N): ";
    string question_2_answer;
    cin>>question_2_answer;//get answer

    // boolean value for saving s parameter data default false
	bool save_s21=false;

    // action take depending on answer
    if((question_2_answer=="n")||(question_2_answer=="no")||
		(question_2_answer=="N")||(question_2_answer=="No")||
		(question_2_answer=="NO")){

        cout<<"\nNot saving.\n";//user answers no, dont save

    }else if((question_2_answer=="y")||(question_2_answer=="yes")||
			(question_2_answer=="Y")||(question_2_answer=="Yes")||
			(question_2_answer=="YES")){

        cout<<"\nsaving data\n";// user answers yes
        save_s21=true;// save data

    }else{
        //unintelligible answer, don't save
		cout<<"\nChoice not understood. Not Saving.\n";
    }
    return save_s21;
}

// save S parameter data to a file
void saveSDataToFile(double *freqs,TwoPortNetwork &tpn,int n_freqs){
    ofstream s_data_file("sdata_albert_patterson_code_sample.txt");
    // headers
    s_data_file<<"Frequency (Hz)\t"
            <<"S11 real\t"
            <<"S11 imag\t"
            <<"S12 real\t"
            <<"S12 imag\t"
            <<"S21 real\t"
            <<"S21 imag\t"
            <<"S22 real\t"
            <<"S22 imag\n";

    // data
    for(int row=0;row<n_freqs;row++){
            s_data_file<<freqs[row]<<"\t";
            for(int col=0; col<8; col++){
                s_data_file<<tpn.get_param_value(row, col)<<"\t";
            }
            s_data_file<<"\n";
    }
}

// save |S21| data to a file
void saveS21DataToFile(double *freqs,double* s21_vals,int n_freqs){
    ofstream s21_data_file("s21data_albert_patterson_code_sample.txt");

    // headers
    s21_data_file<<"Frequency (Hz)\t|S21| (dB)\n";
    //data
    for(int row=0;row<n_freqs;row++){
        s21_data_file<<freqs[row]<<"\t"<<s21_vals[row]<<"\n";
    }
}

// Assess the filter quality for the user
void assessDesign(double insertion_loss, double center_frequency){
    // report on insertion loss and center frequency
    cout<<"\nInsertion loss is "<<insertion_loss<<" dB\nCenter frequency is "
		<<center_frequency/1e9<<" GHz\n";
    // if the design is good:
    if(insertion_loss<=3){
        cout<<"\nPretty good!\n";//low insertion loss
    }else if(insertion_loss>=10){
        cout<<"\nBack to the drawing board!\n";//high insertion loss
    }else{
        cout<<"\nYour design needs work.\n";//moderate insertion loss
    }
}
