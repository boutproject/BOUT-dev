/**
 * @file	BoutIfc.hxx
 *
 * @brief	BOUT++ object implemented through C++
 *
 * @version	$Id: BoutIfc.hxx 1389 2010-06-29 15:27:44Z sveta $
 */

#ifndef BOUT_IFC
#define BOUT_IFC

// configure includes
#ifdef HAVE_CONFIG_H
#include <config.h>
#undef HAVE_CONFIG_H
#endif

// facetsifc
#include <FacetsIfc.h>
#include <Facets0dIfc.h>

// stl includes
#include <string>

/**
 * This basically differentiates a BOUT object from a basic FacetsIfc
 * object.  No new implementation is actually introduced here.
 */
class BoutIfc : public FacetsIfc,
                public Facets0dIfc {

  public:

/**
 * Constructor
 *
 * @param petsc use native uedge petsc?
 * @param mpi use native uedge mpi?
 */
    BoutIfc(int petsc, int mpi);

/**
 * Destructor
 */
    virtual ~BoutIfc();

/**
  Inialize
 */
    virtual int initialize();

/**
 * This is called to set the log file affix
 *
 * @param fname the affix for the filename
 */
    virtual int setLogFile(const std::string& fname);

/**
 * This is called to set the proper communicator
 *
 * @param fint the communicator
 */
    virtual int setMpiComm(long fint);

/**
 * This is called to read additional parameters for uedge
 *
 * @param name the name of the file to read
 */
    virtual int readParams(const std::string& name);

/**
 * This is called to construct all the data required for this solver
 * to function
 */
    virtual int buildData();

/**
 * This is called to construct the solver itself in Uedge
 */
    virtual int buildUpdaters();

/**
 * Get the rank of the interface where it lives
 * @param ifc the name of the interface
 * @param rank on which this interface lives: returned by reference as a parameter
 * @return error code; zero if all is OK
 */
    virtual int getRankOfInterface(const std::string& ifc, size_t& rank) const;

/**
 * Set value of a double parameter
 *
 * @param name name of the variable to set
 * @param value the value to set
 * @return error code; zero if all is OK
 */
    virtual int set0dDouble(const std::string& name, double value);

/**
 * Get value of a double parameter
 * @param name name of the variable to set
 * @param value the returned value
 * @return error code; zero if all is OK
 */
    virtual int get0dDouble(const std::string& name, double& value) const;

/**
 * Set value of integer parameter
 *
 * @param name name of the variable to set
 * @param value the value to set
 * @return error code; zero if all is OK
 */
    virtual int set0dInt(const std::string& name, int value);

/**
 * Get value of integer parameter
 *
 * @param name name of the variable to set
 * @param value the returned value
 * @return error code; zero if all is OK
 */
    virtual int get0dInt(const std::string& name, int& value) const;

/**
 * Set value of string parameter
 *
 * @param name name of the variable to set
 * @param value the value to set
 * @return error code; zero if all is OK
 */
    virtual int set0dString(const std::string& name, const std::string& value);

/**
 * Get value of string parameter
 *
 * @param name name of the variable to set
 * @param value the returned value
 * @return error code; zero if all is OK
 */
    virtual int get0dString(const std::string& name, std::string& value) const;

/**
 * Set a double value at particular location
 *
 * @param name name of the variable
 * @param ndims number of dimensions
 * @param locs locations where to set the value
 * @param val the value to set
 */
    virtual int set0dDoubleAtLocation(const std::string& name, size_t ndims, const double* locs, double val) {return 1;}

/**
 * Get a doublevalue at particular location
 *
 * @param name name of the variable
 * @param ndims number of dimensions
 * @param locs locations where to set the value
 * @param ret the value of the variable
 */
    virtual int get0dDoubleAtLocation(const std::string& name, size_t ndims, const double* locs, double& ret) const {return 1;}


/**
 * Set a int value at particular location
 *
 * @param name name of the variable
 * @param ndims number of dimensions
 * @param locs locations where to set the value
 * @param val the value to set
 */
    virtual int set0dIntAtLocation(const std::string& name, size_t ndims, const double* locs, int val) {return 1;}

/**
 * Get a int value at particular location
 *
 * @param name name of the variable
 * @param ndims number of dimensions
 * @param locs locations where to set the value
 * @param ret the value of the variable
 */
    virtual int get0dIntAtLocation(const std::string& name, size_t ndims, const double* locs, int& ret) const {return 1;}

/**
 * Dump to a file
 * @param file the file name to dump to
 * @return error code; zero if all is OK
 */
    virtual int dumpToFile(const std::string& file) const;

/**
 * Restore from a file
 * @param file the file name to restore from
 * @return error code; zero if all is OK
 */
    virtual int restoreFromFile(const std::string& file);

/**
 * This is called to advance uedge
 *
 * @param t time to which to advance
 */
    virtual int update(double t);

/**
 * Revert the simulation to old time step
 *
 */
    virtual int revert();

/**
 * This is called to complete the uedge object
 */
    virtual int complete();

  protected:
    std::string paramFile;

};

#endif // BOUT_IFC

