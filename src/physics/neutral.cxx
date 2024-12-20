 /**************************************************************************
 * Neutral reaction rates based on Reiter's polynomial fitting formula
 * (http://www.eirene.de/amjuel.pdf)
 *
 **************************************************************************
 * AA is a placeholder for non-hydrogen neutrals (to be improved)
 * by B. Zhu (02/05/2020)
 *
 **************************************************************************/

#include <bout/globals.hxx>
#include <math.h>
#include <bout/neutral.hxx>
#include "bout/output.hxx"

const Field3D iz_rate(const Field3D &tin, BoutReal AA)
{
  Field3D lt, result=0.;

  //output.write("in iz, max and min te=%e,%e.\n", max(tin), min(tin));
  lt = log(tin);
  //output.write("in iz, max and min lt=%e,%e.\n", max(lt), min(lt));

  if (AA == 1) {
    // 2.17 - reaction 2.1.5FJ
    const BoutReal alpha[] = {-0.317385e2, 0.1143818e2, -0.3833998e1, \
           0.7046692, -0.74314862e-1, 0.4153749e-2, -0.9486967e-4};
    for (int n=0; n<7; n++){
      result += alpha[n] * pow(lt,n);
    }
    result = 1.e-6*exp(result); // in [m^3/s]
  }
  else
    output.write("neutral AA=%d hasn't been implemented yet! \n", AA);
  
  //output.write("in iz, max and min result=%e,%e.\n", max(result), min(result));

  return result;
}

const Field3D cx_rate(const Field3D &tin, BoutReal AA)
{
  Field3D lt, result=0.;
  
  //output.write("in cx, max and min te=%e,%e.\n", max(tin), min(tin));
  lt = log(tin);
  //output.write("in cx, max and min lt=%e,%e.\n", max(lt), min(lt));
  
  if (AA == 1) {
    // 3.19 - reaction 3.1.8 with Tn=Ti assumption 
    const BoutReal alpha[] = {-18.3167, 3.793864e-1, -5.1173e-3, 8.5307e-3, \
	  -3.3887e-3, -2.105524e-4, 4.460404e-4, -2.048171e-4, 3.24055e-5, \
	   7.5710e-6, -3.2614e-6, 1.838622e-7, 8.99471e-8, -2.11449e-8, \
      	   2.0666e-9, -9.94415e-011, 1.9350e-12};
    for (int n=0; n<17; n++){
      result += alpha[n] * pow(lt,n);
    }
    result = 1.e-6*exp(result); // in [m^3/s]
  }
  else
    output.write("neutral AA=%d hasn't been implemented yet! \n", AA);

  //output.write("in cx, max and min result=%e,%e.\n", max(result), min(result));

  return result;
}

const Field3D rc_rate(const Field3D &tin, const Field3D &nin, BoutReal AA)
{        
  Field3D lt, ln, result=0.;
  
  //output.write("in rc, max and min te=%e,%e.\n", max(tin), min(tin));
  //output.write("in rc, max and min ni=%e,%e.\n", max(nin), min(nin));
  lt = log(tin);
  ln = log(nin/1.e14);
  //output.write("in rc, max and min lt=%e,%e.\n", max(lt), min(lt));
  //output.write("in rc, max and min ln=%e,%e.\n", max(ln), min(ln));

  if (AA == 1) {
    // 4.4 - reaction 2.1.8JH
    const BoutReal alpha[9][9] = { \
	   {-2.855728479302e1,   3.488563234375e-2, -2.799644392058e-2, \ 
             1.209545317879e-2, -2.436630799820e-3,  2.837893719800e-4, \
            -1.886511169084e-5,  6.752155602894e-7, -1.005893858779e-8}, \
	   {-7.664042607917e-1, -3.583233366133e-3, -7.452514292790e-3, \
             2.709299760454e-3, -7.745129766167e-4,  1.142444698207e-4, \
            -9.382783518064e-6,  3.902800099653e-7, -6.387411585521e-9}, \
	   {-4.930424003280e-3, -3.620245352252e-3,  6.958711963182e-3, \
            -2.139257298118e-3,  4.603883706734e-4, -5.991636837395e-5, \
             4.729262545726e-6, -1.993485395689e-7,  3.352589865190e-9}, \
	   {-5.386830982777e-3, -9.532840484460e-4,  4.631753807534e-4, \
            -5.371179699661e-4,  1.543350502150e-4, -2.257565836876e-5, \
             1.730782954588e-6, -6.618240780594e-8,  1.013364275013e-9}, \
	   {-1.626039237665e-4,  1.888048628708e-4,  1.288577690147e-4, \
            -1.634580516353e-5, -9.601036952725e-6,  3.425262385387e-6, \
            -4.077019941998e-7,  2.042041097083e-8, -3.707977721109e-10},\
	   { 6.080907650243e-6, -1.014890683861e-5, -1.145028889459e-4, \
             5.942193980802e-5, -1.211851723717e-5,  1.118965496365e-6, \
            -4.275321573501e-8,  3.708616111085e-10, 7.068450112690e-12},\
	   { 2.101102051942e-5,  2.245676563601e-5, -2.245624273814e-6, \
            -2.944873763540e-6,  1.002105099354e-6, -1.291320799814e-7, \
             7.786155463269e-9, -2.441127783437e-10, 3.773208484020e-12},\
	   {-2.770717597683e-6, -4.695982369246e-6,  3.250878872873e-6, \
            -9.387290785993e-7,  1.392391630459e-7, -1.139093288575e-8, \
             5.178505597480e-10,-9.452402157390e-12,-4.672724022059e-14},\
	   { 1.038235939800e-7,  2.523166611507e-7, -2.145390398476e-7, \
             7.381435237585e-8, -1.299713684966e-8,  1.265189576423e-9, \
            -6.854203970018e-11, 1.836615031798e-12,-1.640492364811e-14} };
    for (int m=0; m<9; m++){
      for (int n=0; n<9; n++){
        result += alpha[m][n] * pow(lt,m) * pow(ln,n);
      }
    }
    result = 1.e-6*exp(result); // in [m^3/s]
  }
  else
    output.write("neutral AA=%d hasn't been implemented yet! \n", AA);

  //output.write("in rc, max and min result=%e,%e.\n", max(result), min(result));

  return result;
}

const Field3D cx_sect(const Field3D &tin, BoutReal AA)
{
  Field3D lt, result=0.;
  
  lt = log(tin);
  if (AA == 1) {
    // 1.2.3 - reaction 0.1V
    const BoutReal alpha[] = {-3.353420922048e1, -3.522409780724e-1, -3.587214262651e-2, \
          -4.282561006823e-3, -3.230618998917e-4, -4.34317369894e-5, \
          -1.753965583282e-5, -4.580920664987e-7, 3.738689325195e-7};
    for (int n=0; n<9; n++){
      result += alpha[n] * pow(lt,n);
    }
    result = 1.e-4*exp(result); // in [m^2]
  }
  else
    output.write("neutral AA=%d hasn't been implemented yet! \n", AA);

  //output.write("in cx_set, max and min result=%e,%e.\n", max(result), min(result));
  return result;
}

