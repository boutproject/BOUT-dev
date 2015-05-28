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

class RKScheme;

#ifndef __RKSCHEME_H__
#define __RKSCHEME_H__

#include <bout_types.hxx>
#include <options.hxx>
#include <utils.hxx>

#include <string>
#include <vector>
#include <list>

using std::string;
using std::vector;

#define RKSchemeType const char*
#define RKSCHEME_RKF45       "rkf45"
#define RKSCHEME_CASHKARP    "cashkarp"

class RKScheme {
 public:

  RKScheme(Options *opts = NULL); //Options picks the scheme, pretty much everything else is automated
  ~RKScheme();

  //Finish generic initialisation
  void init(const int nlocal);

  //Get the time at given stage
  BoutReal setCurTime(const BoutReal timeIn, const BoutReal dt, const int curStage);

  //Get the state vector at given stage
  void setCurState(const BoutReal *start, BoutReal *out, const int nlocal, 
		   const int curStage, const BoutReal dt);

  //Calculate the two updated states
  void setOutputStates(const BoutReal *start, BoutReal *resultFollow, 
		       BoutReal *resultAlt, const int nlocal, const BoutReal dt);

  //Update the timestep
  virtual BoutReal updateTimestep(const BoutReal dt, const BoutReal rtol, 
				  const BoutReal err);

  //Returns the string name for the given scheme
  virtual string getType(){return label;};

  //Returns the number of stages for the current scheme
  int getStageCount(){return numStages;};

  //The intermediate stages
  BoutReal **steps;

 protected:
  //Information about scheme
  bool followHighOrder; //If true the recommended solution is the higher order one.
  string label;
  int numStages; //Number of stages in the scheme
  int order; //Order of scheme

  //The Butcher Tableau
  BoutReal **stageCoeffs;
  BoutReal **resultCoeffs;
  BoutReal *timeCoeffs;

 private:
};

#endif // __RKSCHEME_H__
