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

#define RKSchemeType const char*
#define RKSCHEME_RKF45       "rkf45"
#define RKSCHEME_CASHKARP    "cashkarp"
#define RKSCHEME_RK4         "rk4"
#define RKSCHEME_RKF34       "rkf34"

class RKScheme {
 public:

  //Options picks the scheme, pretty much everything else is automated
  RKScheme(Options *opts = nullptr);
  virtual ~RKScheme() = default;

  //Finish generic initialisation
  void init(int nlocalIn, int neqIn, bool adaptiveIn, BoutReal atolIn,
            BoutReal rtolIn, Options *options = nullptr);

  //Get the time at given stage
  BoutReal setCurTime(BoutReal timeIn,BoutReal dt,int curStage);

  //Get the state vector at given stage
  virtual void setCurState(const Array<BoutReal> &start, Array<BoutReal> &out,int curStage, 
			   BoutReal dt);

  //Calculate the output state and return the error estimate (if adaptive)
  virtual BoutReal setOutputStates(const Array<BoutReal> &start,BoutReal dt, Array<BoutReal> &resultFollow);

  //Update the timestep
  virtual BoutReal updateTimestep(BoutReal dt,BoutReal err);

  //Returns the string name for the given scheme
  virtual std::string getType(){return label;};

  //Returns the number of stages for the current scheme
  int getStageCount(){return numStages;};

  //Returns the number of orders for the current scheme
  int getNumOrders(){return numOrders;};

  //The intermediate stages
  Matrix<BoutReal> steps;

 protected:
  //Information about scheme
  bool followHighOrder; //If true the recommended solution is the higher order one.
  std::string label;
  int numStages; //Number of stages in the scheme
  int numOrders; //Number of orders in the scheme
  int order; //Order of scheme

  //The Butcher Tableau
  Matrix<BoutReal> stageCoeffs;
  Matrix<BoutReal> resultCoeffs;
  Array<BoutReal> timeCoeffs;
  
  Array<BoutReal> resultAlt;

  int nlocal;
  int neq;
  BoutReal atol;
  BoutReal rtol;
  bool adaptive;

  BoutReal dtfac;

  virtual BoutReal getErr(Array<BoutReal> &solA, Array<BoutReal> &solB);

  virtual void constructOutput(const Array<BoutReal> &start,BoutReal dt, 
			       int index, Array<BoutReal> &sol);

  virtual void constructOutputs(const Array<BoutReal> &start,BoutReal dt, 
				int indexFollow,int indexAlt,
				Array<BoutReal> &solFollow, Array<BoutReal> &solAlt);

 private:
  void verifyCoeffs();
  void printButcherTableau();
  void zeroSteps();
};

#endif // __RKSCHEME_H__
