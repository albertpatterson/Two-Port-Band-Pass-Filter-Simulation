#ifndef TWOPORTNETWORK_H_INCLUDED
#define TWOPORTNETWORK_H_INCLUDED

#include <algorithm>
#include <vector>
#include <fstream>
#include "complexDoubleArray.h"

/*
Abstract class to build and model two port network circuits from resistors,
capacitors and inductors. The simulation can be performed for any frequency
range and with any number of data points using the derived classes.
The s parameters can then be saved.
*/

template<const unsigned int N_FREQS_>
class TwoPortNetwork{
    private:
        const double PI_ = 3.1415926535897;

        // used to assign a linear spacing of frequencies
        // input: start frequency of type double and end frequency of
        //          type double
        // output: none, freqs_ will be spaced linearly from the start
        //          frequency to the end frequency
        void linSpaceFreqs(double start_freq,double end_freq);


        // used to calculate the admittance of a resistor, capacitor, or
        // inductor
        // input: component type: 'r' for resistor, 'c' for capacitor or
        //           'l' for inductor,
        //        component value: value of component in Ohm, Farad or Henry
        // output: complexDoubleArray of admittance values for each frequency
        //          value
        void fillAdmittance(char component_type, double component_val,
                            complexDoubleArray<N_FREQS_>& cda);

        // used to calculate the impedance of a resistor, capacitor, or
        // inductor
        // input: component type: 'r' for resistor, 'c' for capacitor or
        //           'l' for inductor,
        //        component value: value of component in Ohm, Farad or Henry
        // output: complexDoubleArray of impedance values for each frequency
        //          value
        void fillImpedance(char component_type, double component_val,
                           complexDoubleArray<N_FREQS_>& cda);

    protected:

        // frequencies used in simulation in Hz
        // will be assigned by constructor
        array<double,N_FREQS_>* freqs_{new array<double,N_FREQS_>};

        // values of two port network parameters
        // will be assigned by constructor
        array<complexDoubleArray<N_FREQS_>,4> *params_{
                                new array<complexDoubleArray<N_FREQS_>,4>};

    public:

        // construct a passive component for a two-port network
        // input: start_frequency: start frequency of the simulation
        //        end_freq: last frequency of the simulation
        //        component type: 'r','c' or 'l' for resistor, capacitor or
        //          inductor, respectively
        //        component_val: value of component in Ohm, Farad, or Henry
        //        is_series: true for a series component, false for a shunt
        //          component
        TwoPortNetwork(double start_freq, double end_freq, char component_type,
                       double component_val,bool is_series);


        // construct by copying another two port network
        // input: tpn_a: two port network to be copied
        TwoPortNetwork(const TwoPortNetwork<N_FREQS_>& tpn_a){
            *freqs_=*(tpn_a.freqs_);
            *params_=*(tpn_a.params_);
        };

        // destroy two port network
        ~TwoPortNetwork();

        // Create a new two port network by cascading a two port network onto
        //      another two port network
        // input: tpn_a: two port network
        // output: two port network modeling the cascaded connection of this
        //          and tpn_a
        TwoPortNetwork<N_FREQS_> cascade(TwoPortNetwork<N_FREQS_>& tpn_a);


        // cascading a two port network onto another two port network without
        // creating a new twq port network
        // input: tpn_a: two port network
        // output: none, the parameters of this will be altered;
        void cascadeInPlace(TwoPortNetwork<N_FREQS_>& tpn_a);


        // get the s parameters for a two port network
        // input: r0: termination impedance of circuit, default 50 ohm
        // output: an array of 4 complexDoubleArrays for the S11, S12, S21 and
        //          S22, respectively
        array<complexDoubleArray<N_FREQS_>,4> get_s_params(double r0=50);


        // get a single s parameters for a two port network
        // input:  row: row in the table of s parameters
        //          col: column in the table of s parameters
        //          r0: termination impedance of circuit, default 50 ohm
        // output: a complexDouble value for the s parameter
        complexDouble get_s_param(const unsigned int row,
                                  const unsigned int col,double r0);


        //get forward transmission in dB
        // input: r0: termination impedance of circuit, default 50 ohm
        // output: an array of doubles representing the magnitude of S21 in dB
        array<double,N_FREQS_> get_s21_dB(double r0=50);

        // save the s parameters of the two port network
        // input: filename: name for the file
        //        r0: termination impedance, default 50 Ohm
        void saveSParams(string filename, double r0=50);

        // abstract function to find relevant filter properties
        // be implemented in subclasses
        virtual vector<double> get_filter_props(const double r0=50)=0;
};


// used to assign a linear spacing of frequencies
// input: start frequency of type double and end frequency of
//          type double
// output: none, freqs_ will be spaced linearly from the start
//          frequency to the end frequency
template<const unsigned int N_FREQS_>
void TwoPortNetwork<N_FREQS_>::linSpaceFreqs(double start_freq,
                                             double end_freq){

    // spacing between frequencies
    double delt=(end_freq-start_freq)/(N_FREQS_-1);

    int delt_mult=0;

    for(double &freq : *(this->freqs_)){
        freq=start_freq+(delt_mult++)*delt;// increment with delt
    }

}



// used to calculate the admittance of a resistor, capacitor, or
// inductor
// input: component type: 'r' for resistor, 'c' for capacitor or
//           'l' for inductor,
//        component value: value of component in Ohm, Farad or Henry
//        cda: the complexDoubleArray to be filled
// output: none cda will be filled with admittance values
template<const unsigned int N_FREQS_>
void TwoPortNetwork<N_FREQS_>::fillAdmittance(char component_type,
                                            double component_val,
                                            complexDoubleArray<N_FREQS_>& cda){

    // check if component is a capacitor
    if (component_type=='c'){

        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(0,2*PI_*freqs_->at(i)*component_val);
        }

    // check if component is an inductor
    }else if(component_type=='l'){

        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(0,-1/(2*PI_*freqs_->at(i)*component_val));
        }

    // check if component is a resistor
    }else if(component_type=='r'){

        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(1/component_val,0);
        }

    }else{

        // any other value for the component type is not correct
        cout<<"\nerror, incorrect component type: use c for capacitor,"<<
                "l for inductor or r for resistor.Results not valid.\n";
        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(0,0);
        }

    }

}


// used to calculate the impedance of a resistor, capacitor, or
// inductor
// input: component type: 'r' for resistor, 'c' for capacitor or
//           'l' for inductor,
//        component value: value of component in Ohm, Farad or Henry
//        cda: the complexDoubleArray to be filled
// output: none cda will be filled with impedance values
template<const unsigned int N_FREQS_>
void TwoPortNetwork<N_FREQS_>::fillImpedance(char component_type,
                                            double component_val,
                                            complexDoubleArray<N_FREQS_>& cda){

    // check if component is a capacitor
    if (component_type=='c'){

        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(0,-1/(2*PI_*freqs_->at(i)*component_val));
        }

    // check if component is an inductor
    }else if(component_type=='l'){

        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(0,2*PI_*freqs_->at(i)*component_val);
        }

    // check if component is a resistor
    }else if(component_type=='r'){

        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(component_val,0);
        }

    }else{

        // any other value for the component type is not correct
        cout<<"\nerror, incorrect component type: use c for capacitor,"<<
                "l for inductor or r for resistor. Results not valid.\n";
        for(unsigned int i=0;i<N_FREQS_;i++){
            cda[i]=complexDouble(INFINITY,INFINITY);
        }

    }
}



// construct and initialize parameters based from a basic circuit
// element
// input: start_frequency: start frequency of the simulation
//        end_freq: last frequency of the simulation
//        component type: 'r','c' or 'l' for resistor, capacitor or
//          inductor, respectively
//        component_val: value of component in Ohm, Farad, or Henry
//        is_series: true for a series component, false for a shunt
//          component
template<const unsigned int N_FREQS_>
TwoPortNetwork<N_FREQS_>::TwoPortNetwork(double start_freq, double end_freq,
                               char component_type,double component_val,
                               bool is_series){

    // assign freqs_
    linSpaceFreqs(start_freq, end_freq);

    // assign params_
    if (is_series){

            params_->at(0).fill(complexDouble(1,0));
            fillImpedance(component_type,component_val,params_->at(1));
            params_->at(2).fill(complexDouble(0,0));
            params_->at(3).fill(complexDouble(1,0));

    }else{

            params_->at(0).fill(complexDouble(1,0));
            params_->at(1).fill(complexDouble(0,0));
            fillAdmittance(component_type,component_val,params_->at(2));
            params_->at(3).fill(complexDouble(1,0));

    }
}

// Destroy two port network
template<const unsigned int N_FREQS_>
TwoPortNetwork<N_FREQS_>::~TwoPortNetwork(){

    delete freqs_;
    delete[] params_;

}

// cascading a two port network onto another two port network without
// creating a new twq port network
// input: tpn_a: two port network
// output: none, the parameters of this will be altered;
template<const unsigned int N_FREQS_>
void TwoPortNetwork<N_FREQS_>::cascadeInPlace(TwoPortNetwork<N_FREQS_>& tpn_a){

    complexDouble t11, t12, t21, t22;

    // check for frequency mismatch
    if(((this->freqs_)->at(0)==tpn_a.freqs_->at(0))&&
        (((this->freqs_)->at(N_FREQS_-1))==tpn_a.freqs_->at(N_FREQS_-1))){

        // cascade elements
        for(unsigned int i =0; i<N_FREQS_; i++){

            t11=this->params_->at(0)[i]*(tpn_a.params_->at(0)[i])+
                            this->params_->at(1)[i]*tpn_a.params_->at(2)[i];

            t12=this->params_->at(0)[i]*tpn_a.params_->at(1)[i]+
                            this->params_->at(1)[i]*tpn_a.params_->at(3)[i];

            t21=this->params_->at(2)[i]*tpn_a.params_->at(0)[i]+
                            this->params_->at(3)[i]*tpn_a.params_->at(2)[i];

            t22=this->params_->at(2)[i]*tpn_a.params_->at(1)[i]+
                            this->params_->at(3)[i]*tpn_a.params_->at(3)[i];

            this->params_->at(0)[i]=t11;
            this->params_->at(1)[i]=t12;
            this->params_->at(2)[i]=t21;
            this->params_->at(3)[i]=t22;

        }
    }else{

        // if frequencies do not match, cannot cascade elements
        cout<<"\nFrequency Mismatch. No change made.\n";

    }
}

// get a single s parameters for a two port network
// input:  row: row in the table of s parameters
//          col: column in the table of s parameters
//          r0: termination impedance of circuit, default 50 ohm
// output: a complexDouble value for the s parameter
template<const unsigned int N_FREQS_>
complexDouble TwoPortNetwork<N_FREQS_>::get_s_param(const unsigned int row,
                                                    const unsigned int col,
                                                    double r0){

    complexDouble s_param;

    if(col==0){// Calculate S11

        s_param=(params_->at(0)[row]+
                 params_->at(1)[row]/
                    complexDouble(r0,0)-
                 params_->at(2)[row]*
                 complexDouble(r0,0)-
                 params_->at(3)[row])/
                        (params_->at(0)[row]+
                         params_->at(1)[row]/
                            complexDouble(r0,0)+
                         params_->at(2)[row]*
                         complexDouble(r0,0)+
                         params_->at(3)[row]);

    }else if(col==1){//Calculate S12

        s_param=complexDouble(2,0)*
                    (params_->at(0)[row]*
                     params_->at(3)[row]-
                     params_->at(1)[row]*
                     params_->at(2)[row])/
                        (params_->at(0)[row]+
                         params_->at(1)[row]/
                         complexDouble(r0,0)+
                         params_->at(2)[row]*
                         complexDouble(r0,0)+
                         params_->at(3)[row]);

    }else if(col==2){//Calculate S21

        s_param=complexDouble(2,0)/
                    (params_->at(0)[row]+
                     params_->at(1)[row]/
                        complexDouble(r0,0)+
                     params_->at(2)[row]*
                     complexDouble(r0,0)+
                     params_->at(3)[row]);

    }else if (col==3){//Calculate S22

        s_param=(complexDouble(-1,0)*
                 params_->at(0)[row]+
                 params_->at(1)[row]/
                    complexDouble(r0,0)-
                 params_->at(2)[row]*
                 complexDouble(r0,0)+
                 params_->at(3)[row])/
                    (params_->at(0)[row]+
                     params_->at(1)[row]/
                        complexDouble(r0,0)+
                     params_->at(2)[row]*
                     complexDouble(r0,0)+
                     params_->at(3)[row]);

    }else{

        cout<<"\ncol must be 0, 1, 2 or 3 for S11, S12, S21 or S22,"
            <<" respectively.\n";
        s_param=INFINITY;

    }

    return s_param;
}

//get forward transmission in dB
// input: r0: termination impedance of circuit, default 50 ohm
// output: an array of doubles representing the magnitude of S21 in dB
template<const unsigned int N_FREQS_>
array<double,N_FREQS_> TwoPortNetwork<N_FREQS_>::get_s21_dB(double r0){

//   s21=2/(a11+a12/r0+a21*r0+a22)
    complexDoubleArray<N_FREQS_> s_21s=2/(params_->at(0)+params_->at(1)/r0+params_->at(2)*r0+params_->at(3));
    array<double,N_FREQS_> s_21_dbs;
    for(unsigned int i=0; i<N_FREQS_; i++){
        s_21_dbs[i]=20*log10(abs(s_21s[i]));
    }
    return s_21_dbs;
}

// save the s parameters of the two port network
// input: filename: name for the file
//        r0: termination impedance, default 50 Ohm
template<const unsigned int N_FREQS_>
void TwoPortNetwork<N_FREQS_>::saveSParams(string filename,double r0){

    ofstream datafile;
    datafile.open(filename);

    //headers
    datafile<<"Frequencies (Hz)\t";
    datafile<<"S11 Real\t";
    datafile<<"S11 Imag\t";
    datafile<<"S12 Real\t";
    datafile<<"S12 Imag\t";
    datafile<<"S21 Real\t";
    datafile<<"S21 Imag\t";
    datafile<<"S22 Real\t";
    datafile<<"S22 Imag\n";

    // write data
    for(unsigned int row=0; row<N_FREQS_; row++){
            datafile<<freqs_->at(row)<<"\t";
            for(unsigned int col=0; col<4; col++){
                datafile<<get_s_param(row,col,r0).real()<<"\t";
                datafile<<get_s_param(row,col,r0).imag()<<"\t";
            }
            datafile<<"\n";
    }

}

#endif // TWOPORTNETWORK_H_INCLUDED
