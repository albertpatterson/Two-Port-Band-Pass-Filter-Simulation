#ifndef HIGHPASSFILTER_H_INCLUDED
#define HIGHPASSFILTER_H_INCLUDED

#include "TwoPortNetwork.h"

/*
Class to build and model two port Hi-pass filter circuits from resistors,
capacitors and inductors. The simulation can be performed for any frequency
range and with any number of data points. The cutoff frequency can
then be found. Inherits from two port filter class
*/

template<const unsigned int N_FREQS_>
class HighPassFilter : public TwoPortNetwork<N_FREQS_>{

    public:

        // construct by copying another two port network
        // input: tpn: two port network to be copied
        HighPassFilter(const TwoPortNetwork<N_FREQS_>& tpn)
                                :TwoPortNetwork<N_FREQS_>(tpn){};

        // construct a passive component for a hi pass filter
        // input: start_frequency: start frequency of the simulation
        //        end_freq: last frequency of the simulation
        //        component type: 'r','c' or 'l' for resistor, capacitor or
        //          inductor, respectively
        //        component_val: value of component in Ohm, Farad, or Henry
        //        is_series: true for a series component, false for a shunt
        //          component
        HighPassFilter(double start_freq, double end_freq, char component_type,
                           double component_val,bool is_series):
            TwoPortNetwork<N_FREQS_>(start_freq, end_freq, component_type,
                           component_val,is_series){};

        // get the properties of the filter
        // input: termination impedance, default 50 ohm
        // output: cutoff frequency
        vector<double> get_filter_props(double r0=50);
};


// get the properties of the filter
// input: termination impedance, default 50 ohm
// output: cutoff frequency
template<const unsigned int N_FREQS_>
vector<double> HighPassFilter<N_FREQS_>::get_filter_props(double r0){

    vector<double> props;
    double s21db=INFINITY;
    unsigned int i=0;

    // find the first frequency for which |S21| >= -3 dB
    do{
        s21db=20*log10(abs(
                        // calculate s21
                        (complexDouble(2,0)/(this->params_->at(0)[i]
                            +this->params_->at(1)[i]/r0
                            +this->params_->at(2)[i]*r0
                            +this->params_->at(3)[i]))
                        ));
        i++;
    }while((s21db<-3)&&(i<N_FREQS_));

    // add cutoff frequency to props if it was found
    if(s21db<-3){
        cout<<"\nHigh-pass filter cutoff frequency not found.\n";
        props.push_back(INFINITY);
    }else{
        props.push_back(this->freqs_->at(i-1));
    }

    return props;
}

#endif // HIGHPASSFILTER_H_INCLUDED
