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
  
  //Takes current state of variables and returns the evolved vars from both High and Low order schemes
  //along with the embedded error estimate. Note that resultFollow is the result this scheme suggests
  //we keep whilst resultAlt is the other. Some schemes use local extrapolation, meaning we keep the
  //higher order result, whilst others suggest we should keep the lower order one.
  void take_step(BoutReal curtime, BoutReal dt, BoutReal *start, BoutReal *resultFollow, BoutReal *resultAlt, BoutReal errEst);

  //Returns the string name for the given scheme
  virtual string getType(){return label;};

 private:

  void debugPrint(); //Prints the cofficients -- only used in testing

  //Coefficient setup (sets the Butcher Tableau)
  void setStageCoeffs(int stage, const vector<BoutReal> coeffs); //Set the stageCoeffs at a given stage
  void setResultCoeffs(int stage, const BoutReal highOrder, const BoutReal lowOrder); //Set the result coeffs as given stage
  void setTimeCoeffs(int stage, const BoutReal coeff); //Set the time coefficient at given stage

  //Sets which order result we prefer
  void setFollowHighOrder(const bool followHigh = true);

  //The deconstructure Butcher Tableau
  vector< vector<BoutReal> > stageCoeffs; //[numStage][numCoeff] // a_i,j
  vector< vector<BoutReal> > resultCoeffs; //[numStage][order] //b_i
  vector<BoutReal> timeCoeffs; //[numStage] //c_j

  //The intermediate stages
  vector< vector<BoutReal> > steps; //[numstage][nlocal]

  //Information about scheme
  int numStages; //Number of stages in the scheme
  bool followHighOrder; //If true the recommended solution is the higher order one.
  string type; //What type of scheme am I?
  string label;
};

#endif // __RKSCHEME_H__
