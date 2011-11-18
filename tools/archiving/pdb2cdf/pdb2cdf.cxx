/*******************************************************
 * PDB2CDF
 * 
 * Convert PDB files to netCDF
 *******************************************************/

#include <stdio.h>
#include <string.h>

// The PDB library (C)
#include "pdb.h"

// netCDF C++ library
#include <netcdfcpp.h>

// Dimension
struct TDim {
  char* label; // Label for the dimension
  int size;    
  int min, max; // Minimum, maximum index
  NcDim* nDim; // netCDF dimension
};

// List of dimensions. Handles up to 3
static TDim dimlist3d[] = {{"x", 0},
			   {"y", 0},
			   {"z", 0}};

// Special case for 4D
static TDim dimlist4d[] = {{"t", 0},
			   {"x", 0},
			   {"y", 0},
			   {"z", 0}};

int main(int argc, char** argv)
{
  TDim *dimlist;
  
  if(argc < 2) {
    fprintf(stderr, "Useage: %s file1 file2 ...\n", argv[0]);
    return 1;
  }

  for(int i=1;i<argc;i++) { // Go through each PDB file
    char *inname = argv[i];
    PDBfile *in;
    
    // Open input file

    if((in = PD_open(inname, "r")) == NULL) {
      fprintf(stderr, "ERROR: Could not open input file '%s'\n", inname);
      return 1;
    }

    // Get list of variables
    
    int nvars;
    char **var_names = PD_ls(in, NULL, NULL, &nvars);
    if((var_names == (char**) NULL) || (nvars < 1)) {
      fprintf(stderr, "ERROR: No variables\n");
      return 1;
    }
    
    // Create a filename for the output
    char *outname;
    int len = strlen(inname);
    if(len < 5) {
      // Not long enough - just append .cdl
      outname = new char[len+5];
      snprintf(outname, len+5, "%s%s", inname, ".cdl");
    }else {
      // Replace '.pdb' by '.cdl'
      outname = new char[len+1];
      strncpy(outname, inname, len+1);
      strncpy(outname+len-3, "cdl", 3);
    }

    // Open output file, overwriting if it exists
    NcFile dataFile(outname, NcFile::Replace);

    if (!dataFile.is_valid()) {
      fprintf(stderr, "ERROR: Could not open output file '%s'\n", outname);
      return 1;
    }

    printf("Converting %s -> %s\n", inname, outname);

    // Find the variable with the largest number of dimensions
    printf("\tAnalysing dimensions...");
    int maxdims = -1;  // maximum number of variables
    int maxdimvar = 0; // variable with the most dimensions
    
    syment *ep;  // PDB query types
    dimdes* dims;

    char *varname;
    for(int var=0;var<nvars;var++) {
      varname = var_names[var];
      // Query size of the variable
      
      if((ep = PD_query_entry(in, varname, NULL)) == (syment*) NULL) {
	fprintf(stderr, "Error querying variable %s\n", varname);
	return 1;
      }
      dims = PD_entry_dimensions(ep);
      int nd = 0; /* Count number of dimensions */
      while(dims != (dimdes*) NULL) {
	nd++;
	dims = dims->next;
      }
      if(nd > maxdims) {
	maxdims = nd;
	maxdimvar = var;
      }
    }

    printf("%d dimensions\n", maxdims);
    
    if(maxdims < 4) {
      dimlist = dimlist3d;
    }else
      dimlist = dimlist4d;
    
    if(maxdims > 4) {
      fprintf(stderr, "ERROR: Can only handle up to 4D variables\n");
      return 1;
    }
    
    // Get the size of each dimension
    varname = var_names[maxdimvar];
    if((ep = PD_query_entry(in, varname, NULL)) == (syment*) NULL) {
      fprintf(stderr, "Error querying variable %s\n", varname);
      return 1;
    }
    dims = PD_entry_dimensions(ep);
    
    for(int d=0;d<maxdims;d++) {
      dimlist[d].min = dims->index_min;
      dimlist[d].max = dims->index_max;
      dimlist[d].size = dims->index_max - dims->index_min + 1;
      
      // Create a netCDF dimension
      dimlist[d].nDim = dataFile.add_dim(dimlist[d].label, dimlist[d].size);
      
      printf("\t\t%s: %d -> %d (%d)\n", dimlist[d].label, dimlist[d].min, dimlist[d].max, dimlist[d].size);

      dims = dims->next;
    }

    // Go through each variable
    for(int var=0;var<nvars;var++) {
      varname = var_names[var];
      if((ep = PD_query_entry(in, varname, NULL)) == (syment*) NULL) {
	fprintf(stderr, "Error querying variable %s\n", varname);
	return 1;
      }
      
      printf("\t%s", varname);

      // Get dimensions
      int nd = 0; // Number of dimensions
      int vardim[4]; // Indices in the dimlist array
      int lastdim = -1; // Always assume indices have the same order
      int varsize = 1; // Number of elements

      bool gotdims = true;

      dims = PD_entry_dimensions(ep);
      while(dims != (dimdes*) NULL) {
	int min, max;
	min = dims->index_min;
	max = dims->index_max;
	
	varsize *= max - min + 1;
	
	vardim[nd] = -1;
	// Find which dimension this corresponds to
	for(int d = lastdim+1;d<maxdims;d++) {
	  if((dimlist[d].min == min) && (dimlist[d].max == max)) {
	    if(lastdim == -1) {
	      printf("[%s", dimlist[d].label);
	    }else
	      printf(",%s", dimlist[d].label);
	    vardim[nd] = d;
	    lastdim = d;
	    break;
	  }
	}
	if(vardim[nd] == -1) {
	  // Not an existing dimension. Should create a new dimension
	  fprintf(stderr, "ERROR: %s has an unrecognised %d dimension [%d -> %d]\n", varname, nd+1, min, max);
	  gotdims = false;
	  break;
	}
	
	nd++;
	dims = dims->next;
      }
      
      if(!gotdims)
	continue; // Skip this variable

      if(lastdim != -1)
	printf("] (%d elements) ", varsize);

      // Now know number of dimensions nd, and a list of dimension indices vardim

      // Get variable type
      char *typ;
      typ = PD_entry_type(ep);

      printf(" Type: %s ", typ);
      fflush(stdout);

      if(strcasecmp(typ, "integer") == 0) {
	
	int *idata = new int[varsize];
	
	// Read the data from the PDB file
	if (PD_read_as(in, varname, "integer", idata) == FALSE) {
	  fprintf(stderr, "\t\tWARNING: Could not read variable. Ignoring\n");
	  continue;
	}
	switch(nd) {
	case 0: {
	  // Add a 0-D variable to the netCDF file
	  NcVar *ncdata = dataFile.add_var(varname, ncInt);
	  
	  // Write data
	  ncdata->put(idata);
	  break;
	}
	case 1: {
	  NcVar *ncdata = dataFile.add_var(varname, ncInt, dimlist[vardim[0]].nDim);
	  ncdata->put(idata, varsize);
	  break;
	}
	case 2: {
	  NcVar *ncdata = dataFile.add_var(varname, ncInt, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim);
	  ncdata->put(idata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size);
	  break;
	}
	case 3: {
	  NcVar *ncdata = dataFile.add_var(varname, ncInt, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim,
					   dimlist[vardim[2]].nDim);
	  ncdata->put(idata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size,
		      dimlist[vardim[2]].size);
	  break;
	}
	case 4: {
	  NcVar *ncdata = dataFile.add_var(varname, ncInt, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim,
					   dimlist[vardim[2]].nDim,
					   dimlist[vardim[3]].nDim);
	  ncdata->put(idata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size,
		      dimlist[vardim[2]].size,
		      dimlist[vardim[3]].size);
	  break;
	}
	default: {
	  fprintf(stderr, "\t\tWARNING: Cannot yet handle %d-D integers. Ignoring\n", nd); 
	}
	}
	delete[] idata;

	//////////////////////////////////////////////////////////////////
	
      }else if(strcasecmp(typ, "float") == 0) {
	
	float *fdata = new float[varsize];
	
	// Read the data from the PDB file
	if (PD_read_as(in, varname, "float", fdata) == FALSE) {
	  fprintf(stderr, "\t\tWARNING: Could not read variable. Ignoring\n");
	  continue;
	}
	switch(nd) {
	case 0: {
	  NcVar *ncdata = dataFile.add_var(varname, ncFloat);
	  
	  // Write data
	  ncdata->put(fdata);
	  break;
	}
	case 1: {
	  NcVar *ncdata = dataFile.add_var(varname, ncFloat, dimlist[vardim[0]].nDim);
	  ncdata->put(fdata, varsize);
	  break;
	}
	case 2: {
	  NcVar *ncdata = dataFile.add_var(varname, ncFloat, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim);
	  ncdata->put(fdata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size);
	  break;
	}
	case 3: {
	  NcVar *ncdata = dataFile.add_var(varname, ncFloat, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim,
					   dimlist[vardim[2]].nDim);
	  ncdata->put(fdata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size,
		      dimlist[vardim[2]].size);
	  break;
	}
	case 4: {
	  NcVar *ncdata = dataFile.add_var(varname, ncFloat, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim,
					   dimlist[vardim[2]].nDim,
					   dimlist[vardim[3]].nDim);
	  ncdata->put(fdata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size,
		      dimlist[vardim[2]].size,
		      dimlist[vardim[3]].size);
	  break;
	}
	default: {
	  fprintf(stderr, "\t\tWARNING: Cannot yet handle %d-D floats. Ignoring\n", nd); 
	}
	}
	delete[] fdata;
	
	//////////////////////////////////////////////////////////////////
	
      }else if(strcasecmp(typ, "double") == 0) {

	double *ddata = new double[varsize];
	
	// Read the data from the PDB file
	if (PD_read_as(in, varname, "double", ddata) == FALSE) {
	  fprintf(stderr, "\t\tWARNING: Could not read variable. Ignoring\n");
	  continue;
	}
	switch(nd) {
	case 0: {
	  NcVar *ncdata = dataFile.add_var(varname, ncDouble);
	  
	  // Write data
	  ncdata->put(ddata);
	  break;
	}
	case 1: {
	  NcVar *ncdata = dataFile.add_var(varname, ncDouble, dimlist[vardim[0]].nDim);
	  ncdata->put(ddata, varsize);
	  break;
	}
	case 2: {
	  NcVar *ncdata = dataFile.add_var(varname, ncDouble, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim);
	  ncdata->put(ddata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size);
	  break;
	}
	case 3: {
	  NcVar *ncdata = dataFile.add_var(varname, ncDouble, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim,
					   dimlist[vardim[2]].nDim);
	  ncdata->put(ddata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size,
		      dimlist[vardim[2]].size);
	  break;
	}
	case 4: {
	  NcVar *ncdata = dataFile.add_var(varname, ncDouble, 
					   dimlist[vardim[0]].nDim, 
					   dimlist[vardim[1]].nDim,
					   dimlist[vardim[2]].nDim,
					   dimlist[vardim[3]].nDim);
	  ncdata->put(ddata, 
		      dimlist[vardim[0]].size,
		      dimlist[vardim[1]].size,
		      dimlist[vardim[2]].size,
		      dimlist[vardim[3]].size);
	  break;
	}
	default: {
	  fprintf(stderr, "\t\tWARNING: Cannot yet handle %d-D doubles. Ignoring\n", nd); 
	}
	}
	delete[] ddata;
      }else {
	fprintf(stderr, "\tWARNING: '%s' has unrecognised type '%s'. Ignoring\n", varname, typ);
      }
      printf("\n");
    }
    
    delete[] outname; 

    dataFile.close(); // Close the output file. Probably optional.
  }
  
  return 0;
}

