#include <iostream>

#include "BandPassFilter.h"

using namespace std;


/*
Example of a simple simulation of the frequency response of a typical 2-port
AlN, contour-mode resonator and third order filter operating at approximately
1 GHz
*/

int main()
{
    // frequency range to simulate
    const unsigned int n_freqs=10000;// number of frequency values in simulation
    double start_freq=.9e9;//900 MHz initial frequency of simulation
    double end_freq=1.1e9;//1100 MHz final frequency of simulation

    // motional resistance of a resonator
    double motional_resistance_val=20;// 20 ohms
    //motional capacitance of a resonator
	double motional_capacitance_val=4.26e-15;//4.26 fF;
    // motional inductance of a resonator
	double motional_inductance_val=5.95e-6;//5.95 uH;
     //static capacitance of a resonator
	double static_capacitance_val=525e-15;// 525 fF;

    // termination impedance
	double r0=300;

    // create a two port network for the shunt static capacitor
	BandPassFilter<n_freqs> static_capacitor(start_freq,end_freq,'c',
                                            static_capacitance_val,
                                            false);

    // create a two port network for the series motional resistor
	BandPassFilter<n_freqs> motional_resistor(start_freq,end_freq,'r',
                                            motional_resistance_val,
                                            true);

    // create a two port network for the series motional inductor
	BandPassFilter<n_freqs> motional_inductor(start_freq,end_freq, 'l',
                                            motional_inductance_val,
                                            true);

    // create a two port network for the series motional capacitor
	BandPassFilter<n_freqs> motional_capacitor(start_freq,end_freq,'c',
                                                motional_capacitance_val,true);


    // create a resonator by cascading the static capacitor, motional resistor
    // motional inductor, motional capacitor, and the static capacitor again
    BandPassFilter<n_freqs> resonator(static_capacitor);
    resonator.cascadeInPlace(motional_resistor);
    resonator.cascadeInPlace(motional_inductor);
    resonator.cascadeInPlace(motional_capacitor);
    resonator.cascadeInPlace(static_capacitor);


    // create a third order filter by cascading three resonators
    BandPassFilter<n_freqs> third_order_filter(resonator);
        third_order_filter.cascadeInPlace(resonator);
        third_order_filter.cascadeInPlace(resonator);

    // save the s parameters of the resonator, matched to the termination
    // impedance, uncomment below to save
    //resonator.saveSParams("resonator_sps.txt",r0);

    // save the s parameters of the third-order filter,matched to the
    // termination impedance, uncomment below to save
    //third_order_filter.saveSParams("third_order_filter_sps.txt",r0);

    // find the insertion loss and center frequency of the third-order filter
    vector<double> il_f0=third_order_filter.get_filter_props(r0);

    // report results
    cout<<"insertion loss:\t\t"<<il_f0[0]<<" dB"<<"\ncenter frequency:\t"
        <<il_f0[1]/1e9<<"GHz \n\n";

    return 0;
}
