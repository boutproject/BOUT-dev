/**************************************************************************
 * Base class for rk-generic schemes
 *
 **************************************************************************
 * Written by D Dickinson, 2015
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

//////////////////////////////////////////
//// COMMENTS
/*   Originally designed to deal with embedded schemes.
     Doesn't currently support FSAL optimisations.
     Would be nice to replace the coeff arrays with stl containers.
     Could perhaps add a flag to enable "local extrapolation"
*/
//////////////////////////////////////////

class RKScheme;

#ifndef __RKSCHEME_H__
#define __RKSCHEME_H__

#include <bout_types.hxx>
#include <options.hxx>
#include <utils.hxx>

#include <iomanip>
#include <string>
using std::string;
using std::setw;

#define RKSchemeType const char*
#define RKSCHEME_RKF45       "rkf45"
#define RKSCHEME_CASHKARP    "cashkarp"

class RKScheme {
 public:

  RKScheme(Options *opts = NULL); //Options picks the scheme, pretty much everything else is automated
  ~RKScheme();

  //Finish generic initialisation
  void init(const int nlocalIn, const int neqIn, 
	    const BoutReal atolIn, const BoutReal rtolIn);

  //Get the time at given stage
  BoutReal setCurTime(const BoutReal timeIn, const BoutReal dt, const int curStage);

  //Get the state vector at given stage
  void setCurState(const BoutReal *start, BoutReal *out, const int curStage, const BoutReal dt);

  //Calculate the two updated states
  void setOutputStates(const BoutReal *start, BoutReal *resultFollow, 
		       BoutReal *resultAlt, const BoutReal dt);

  //Update the timestep
  virtual BoutReal updateTimestep(const BoutReal dt, const BoutReal err);

  //Returns the string name for the given scheme
  virtual string getType(){return label;};

  //Returns the number of stages for the current scheme
  int getStageCount(){return numStages;};

  //Returns the number of orders for the current scheme
  int getNumOrders(){return numOrders;};

  //The intermediate stages
  BoutReal **steps;

 protected:
  //Information about scheme
  bool followHighOrder; //If true the recommended solution is the higher order one.
  string label;
  int numStages; //Number of stages in the scheme
  int numOrders; //Number of orders in the scheme
  int order; //Order of scheme

  //The Butcher Tableau
  BoutReal **stageCoeffs;
  BoutReal **resultCoeffs;
  BoutReal *timeCoeffs;

 private:
  int nlocal;
  int neq;
  BoutReal atol;
  BoutReal rtol;

  void verifyCoeffs();
  void printButcherTableau();
  void zeroSteps();
};

#endif // __RKSCHEME_H__
