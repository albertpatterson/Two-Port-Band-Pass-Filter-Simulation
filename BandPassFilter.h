#ifndef BANDPASSFILTER_H_INCLUDED
#define BANDPASSFILTER_H_INCLUDED

#include "TwoPortNetwork.h"

/*
Class to build and model two port band-pass filter circuits from resistors,
capacitors and inductors. The simulation can be performed for any frequency
range and with any number of data points. The insertion loss and the center
frequency can then be found. Inherits from two port filter class
*/

template<const unsigned int N_FREQS_>
class BandPassFilter : public TwoPortNetwork<N_FREQS_>{

    public:

        // construct by copying another two port network
        // input: tpn: two port network to be copied
        BandPassFilter(const TwoPortNetwork<N_FREQS_>& tpn)
                                :TwoPortNetwork<N_FREQS_>(tpn){};

        // construct a passive component for a band-pass filter
        // input: start_frequency: start frequency of the simulation
        //        end_freq: last frequency of the simulation
        //        component type: 'r','c' or 'l' for resistor, capacitor or
        //          inductor, respectively
        //        component_val: value of component in Ohm, Farad, or Henry
        //        is_series: true for a series component, false for a shunt
        //          component
        BandPassFilter(double start_freq, double end_freq, char component_type,
                           double component_val,bool is_series):
            TwoPortNetwork<N_FREQS_>(start_freq, end_freq, component_type,
                           component_val,is_series){};

        // get the properties of the filter
        // input: termination impedance, default 50 ohm
        // output: insertion loss, center frequency
        vector<double> get_filter_props(double r0=50);
};

// get the properties of the filter
// input: termination impedance, default 50 ohm
// output: insertion loss, center frequency
template<const unsigned int N_FREQS_>
vector<double> BandPassFilter<N_FREQS_>::get_filter_props(double r0){

    vector<double> props;

    // calculate |S21| in dB
    array<double,N_FREQS_> *s21= new array<double,N_FREQS_>;
    for(unsigned int i=0; i<N_FREQS_; i++){
        s21->at(i)=20*log10(abs(
                            // calculate s21
                            (complexDouble(2,0)/(this->params_->at(0)[i]
                                +this->params_->at(1)[i]/r0
                                +this->params_->at(2)[i]*r0
                                +this->params_->at(3)[i]))
                            ));
    }

    // Find the maximum value nf |S21|
    auto maxval=max_element(s21->begin(),s21->end());

    // add insertion loss to props
    props.push_back(-1*(*maxval));
    // add center frequency to props
    props.push_back((this->freqs_)->at(distance(s21->begin(),maxval)));

    delete s21;
    return props;
}

#endif // BANDPASSFILTER_H_INCLUDED
