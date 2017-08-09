/*
 * mathematica.h
 *
 *  Created on: Dec 11, 2013
 *      Author: yolen
 */

#ifndef MATHEMATICA_H_
#define MATHEMATICA_H_

#include <bout.hxx>
#include <math.h>

BoutReal Power(BoutReal x,BoutReal n) {return pow((double)x,(double)n);}
BoutReal Sin(BoutReal x) {return sin((double)x);}
BoutReal Cos(BoutReal x) {return cos((double)x);}

const BoutReal E = exp(1);


#endif /* MATHEMATICA_H_ */

