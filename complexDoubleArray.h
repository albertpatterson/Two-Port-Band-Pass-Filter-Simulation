#ifndef COMPONENTANALYSIS_H_INCLUDED
#define COMPONENTANALYSIS_H_INCLUDED

#include <complex>
#include <array>

using namespace std;

/*
library to make operations with arrays of complex numbers very convenient
*/

// easy type to use in place of complex<double>
using complexDouble=complex<double>;

// easy type to use in place of array<complex<double>,array_length>
template<const unsigned int array_length> using complexDoubleArray=
                        array<complexDouble,array_length>;

// apply an element by element operation using a complexDoubleArray and a
// complexDouble
// input: complexDoubleArray,
//      complexDouble,
//      operation to apply ex: plus<complexDouble>()
template<const unsigned int array_length,typename operation>
complexDoubleArray<array_length> operate_element_by_element(
                            const complexDoubleArray<array_length> &cd_array1,
                            const complexDoubleArray<array_length> &cd_array2,
                            operation op){

    complexDoubleArray<array_length> return_cd_array;
    // apply operation element by element
    for(unsigned int i=0; i<array_length;i++){
        return_cd_array[i]=op(cd_array1[i],cd_array2[i]);
    }
    return return_cd_array;
}

// apply an element by element operation using a complexDoubleArray and a
// complexDouble
// input: complexDoubleArray,
//      complexDouble,
//      operation to apply ex: plus<complexDouble>()
template<const unsigned int array_length,typename operation>
complexDoubleArray<array_length> operate_element_by_element(
                            const complexDoubleArray<array_length> &cd_array,
                            complexDouble cd, operation op){

    complexDoubleArray<array_length> return_cd_array;
    // apply operation element by element
    for(unsigned int i=0; i<array_length;i++){
        return_cd_array[i]=op(cd_array[i],cd);
    }
    return return_cd_array;
}

// apply an element by element operation using a complexDoubleArray and a
// complexDouble
// input: complexDouble,
//      complexDoubleArray,
//      operation to apply ex: plus<complexDouble>()
template<const unsigned int array_length,typename operation>
complexDoubleArray<array_length> operate_element_by_element(complexDouble cd,
                             const complexDoubleArray<array_length> &cd_array,
                             const operation op){

    complexDoubleArray<array_length> return_cd_array;
    // apply operation element by element
    for(unsigned int i=0; i<array_length;i++){
        return_cd_array[i]=op(cd,cd_array[i]);
    }
    return return_cd_array;
}


//******************************Addition***************************************

// add a complex value to an array of complex values element by element
// input: complex values of type complexDoubleArray<size_t>,
//      complex value of type std::complexDouble
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator +(
                            const complexDoubleArray<array_length> &cd_array1,
                            const complexDoubleArray<array_length> &cd_array2){

    return operate_element_by_element(cd_array1,cd_array2,
                                      plus<complexDouble>());
}


// add a complex value to an array of complex values element by element
// input: complex values of type complexDoubleArray<size_t>,
//      complex value of type std::complexDouble
// output: array of complex values of type
//            complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator +(
                            const complexDoubleArray<array_length> &cd_array,
                            const complexDouble cd){

    return operate_element_by_element(cd_array,cd,plus<complexDouble>());
}

// add a complex value to an array of complex values element by element
// input: complex value of type std::complexDouble,
//      array of complex values of type complexDoubleArray<const unsigned int>
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator +(
                            const complexDouble cd,
                            const complexDoubleArray<array_length> &cd_array){
    return cd_array+cd;
}

// add a real value to an array of complex values element by element
// input: complex values of type complexDoubleArray<const unsigned int>,
//      real value of type double
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator +(
                            const complexDoubleArray<array_length> &cd_array,
                            const double real_value){

        return cd_array+complexDouble(real_value,0);
}

// add a real value to an array of complex values element by element
// input: real value of type double,
//      complex values of type complexDoubleArray<const unsigned int>
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator +(
                            const double real_value,
                            complexDoubleArray<array_length> &cd_array){

        return cd_array+complexDouble(real_value,0);
}

//******************************Addition***************************************


//******************************Subtraction*************************************

// add a complex value to an array of complex values element by element
// input: complex values of type complexDoubleArray<size_t>,
//      complex value of type std::complexDouble
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator -(
                            const complexDoubleArray<array_length> &cd_array1,
                            const complexDoubleArray<array_length> &cd_array2){

    return operate_element_by_element(cd_array1,cd_array2,
                                      minus<complexDouble>());
}

// subtract a complex value from an array of complex values element by element
// input: array of complex values of type
//              complexDoubleArray<const unsigned int>,
//      complex value of type std::complexDouble
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator -(
                            const complexDoubleArray<array_length> &cd_array,
                            const complexDouble cd){
// return element by element operation
    return operate_element_by_element(cd_array,cd,minus<complexDouble>());
}

// subtract an array of complex value from a complex value element by element
// input: complex value of type std::complexDouble,
//      array of complex values of type complexDoubleArray<const unsigned int>
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator -(
                            const complexDouble cd,
                            const complexDoubleArray<array_length> &cd_array){
// return element by element operation
    return operate_element_by_element(cd,cd_array,minus<complexDouble>());
}

// subtract a real value from an array of complex values element by element
// input: array of complex values of type
//          complexDoubleArray<const unsigned int>,
//      real value of type double
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator -(
                            const complexDoubleArray<array_length> &cd_array,
                            const double real_value){

    return cd_array-complexDouble(0,real_value);
}

// subtract an array of complex values from a real value element by element
// input: real value of type double,
//      array of complex values of type complexDoubleArray<const unsigned int>
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator -(
                            const double real_value,
                            const complexDoubleArray<array_length> &cd_array){

    return complexDouble(0,real_value)-cd_array;
}

//******************************Subtraction************************************


//******************************Multiplication*********************************

// add a complex value to an array of complex values element by element
// input: complex values of type complexDoubleArray<size_t>,
//      complex value of type std::complexDouble
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator *(
                            const complexDoubleArray<array_length> &cd_array1,
                            const complexDoubleArray<array_length> &cd_array2){

    return operate_element_by_element(cd_array1,cd_array2,
                                      multiplies<complexDouble>());
}

// multiply an array of complex values by a complex value element by element
// input: complex values of type complexDoubleArray<const unsigned int>,
//      array of complex values of type std::complexDouble
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator *(
                            const complexDoubleArray<array_length> &cd_array,
                            const complexDouble cd){
// return element by element operation
    return operate_element_by_element(cd_array,cd,multiplies<complexDouble>());
}

// multiply a complex value by an array of complex values element by element
// input: complex value of type std::complexDouble,
//      array of complex values of type complexDoubleArray<const unsigned int>
// output: array of complex values of type
//              complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator *(
                            const complexDouble cd,
                            const complexDoubleArray<array_length> &cd_array){

    return cd_array*cd;
}

// multiply a real value by an array of complex values element by element
// input: real of type double,
//      array of complex values of type complexDoubleArray<const unsigned int>
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator *(
                            const double real_value,
                            const complexDoubleArray<array_length> &cd_array){

    return cd_array*complexDouble(real_value,0);
}

// multiply an array of complex values by a real value element by element
// input: array of complex values of type
//              complexDoubleArray<const unsigned int>,
//      real of type double
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator *(
                            const complexDoubleArray<array_length> &cd_array,
                            const double real_value){

    return cd_array*complexDouble(real_value,0);
}

//******************************Multiplication*********************************


//********************************Division*************************************

// add a complex value to an array of complex values element by element
// input: complex values of type complexDoubleArray<size_t>,
//      complex value of type std::complexDouble
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator /(
                            const complexDoubleArray<array_length> &cd_array1,
                            const complexDoubleArray<array_length> &cd_array2){

    return operate_element_by_element(cd_array1,cd_array2,
                                      divides<complexDouble>());
}

// divide an array of complex values by a complex value element by element
// input: array of complex values of type
//          complexDoubleArray<const unsigned int>,
//      complex value of type std::complexDouble
// output: array of complex values of type
//              complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator /(
                            const complexDoubleArray<array_length> &cd_array,
                            const complexDouble cd){
// return element by element operation
    return operate_element_by_element(cd_array,cd,divides<complexDouble>());
}

// divide a complex value by an array of complex values element by element
// input: array of complex value of type std::complexDouble,
//      array of complex values of type complexDoubleArray<const unsigned int>
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator /(
                            const complexDouble cd,
                            const complexDoubleArray<array_length> &cd_array){
    // return element by element operation
    return operate_element_by_element(cd,cd_array,divides<complexDouble>());
}

// divide an array of complex values by a real value element by element
// input: array of complex values of type
//          complexDoubleArray<const unsigned int>,
//      real value of type double
// output: array of complex values of type
//          complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator /(
                            const complexDoubleArray<array_length> &cd_array,
                            const double real_value){

    return cd_array/complexDouble(real_value,0);
}

// divide a real value by an array of complex values element by element
// input: real value of type double,
//      array of complex values of type
//      complexDoubleArray<const unsigned int>
// output: array of complex values of type
//      `complexDoubleArray<const unsigned int>
template<const unsigned int array_length>
complexDoubleArray<array_length> operator /(
                            const double real_value,
                            const complexDoubleArray<array_length> &cd_array){

    return complexDouble(real_value,0)/cd_array;
}

//********************************Division*************************************

#endif // COMPONENTANALYSIS_H_INCLUDED
