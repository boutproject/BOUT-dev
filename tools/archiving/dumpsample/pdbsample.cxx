/*******************************************************
 * PDB Sample
 *
 *
 * Ben Dudson, University of York, July 2009
 *******************************************************/

#include <stdio.h>
#include <string.h>

// The PDB library (C)
#include "pdb.h"

int main(int argc, char** argv) {
  if (argc < 4) {
    fprintf(stderr, "Useage: %s <input file> <output file> <t stride>\n", argv[0]);
    return 1;
  }

  int tstride;
  if (sscanf(argv[3], "%d", &tstride) != 1) {
    fprintf(stderr, "\tERROR: t stride must be an integer\n");
    return 1;
  }

  // Open input file
  PDBfile* in;
  if ((in = PD_open(argv[1], "r")) == NULL) {
    fprintf(stderr, "\tERROR: Could not open input file '%s'\n", argv[1]);
    return 1;
  }

  // Open output file
  PDBfile* out;
  if ((out = PD_open(argv[2], "w")) == NULL) {
    fprintf(stderr, "\tERROR: Could not open output file '%s'\n", argv[2]);
    return 1;
  }

  // Get list of variables

  int nvars;
  char** var_names = PD_ls(in, NULL, NULL, &nvars);
  if ((var_names == (char**)NULL) || (nvars < 1)) {
    fprintf(stderr, "\tERROR: No variables\n");
    return 1;
  }

  // Go through the variables
  char* varname;
  for (int var = 0; var < nvars; var++) {
    varname = var_names[var];

    syment* ep; // PDB query types
    dimdes* dims;

    // Query size of the variable
    if ((ep = PD_query_entry(in, varname, NULL)) == (syment*)NULL) {
      fprintf(stderr, "\tError querying variable %s\n", varname);
      return 1;
    }
    dims = PD_entry_dimensions(ep);
    int nd = 0;      // Count number of dimensions
    int varsize = 1; // Number of elements
    long inds[12];
    while (dims != (dimdes*)NULL) {
      long min, max;
      min = dims->index_min;
      max = dims->index_max;

      if (nd > 3) {
        fprintf(stderr, "\tERROR: Can't handle variable '%s': more than 4D\n", varname);
        return 2;
      }

      inds[3 * nd] = min;
      inds[3 * nd + 1] = max;
      inds[3 * nd + 2] = 1L;

      varsize *= max - min + 1;

      nd++;
      dims = dims->next;
    }

    // Get variable type
    char* typ;
    typ = PD_entry_type(ep);

    if ((strcmp(varname, "t_array") == 0) && (nd == 1)) {
      float* fdata = new float[varsize];

      // Read the data from the PDB file

      inds[2] = tstride;
      int nread;
      if ((nread = PD_read_as_alt(in, varname, "float", fdata, inds)) == 0) {
        fprintf(stderr, "\tWARNING: Could not read t_array. Ignoring\n");
        delete[] fdata;
        continue;
      }

      inds[0] = 0L;
      inds[1] = nread - 1;
      inds[2] = 1L;

      if (PD_write_alt(out, varname, "float", fdata, 1, inds) == FALSE) {
        fprintf(stderr, "\tWARNING: Could not write '%s'. Ignoring\n", varname);
      }

      delete[] fdata;

    } else if (nd == 4) {
      // Reducing time resolution

      if (strcasecmp(typ, "integer") == 0) {
        int* idata = new int[varsize];

        // Read the data from the PDB file
        inds[2] = tstride;
        int nread;
        if ((nread = PD_read_as_alt(in, varname, "integer", idata, inds)) == 0) {
          fprintf(stderr, "\tWARNING: Could not read '%s'. Ignoring\n", varname);
          delete[] idata;
          continue;
        }

        if (nread % (varsize / (inds[1] - inds[0] + 1)) != 0) {
          fprintf(stderr, "ERROR: Accounting error: (%ld, %ld), %d, %d\n", inds[0],
                  inds[1], varsize, nread);
          delete[] idata;
          continue;
        }

        nread = nread / (varsize / (inds[1] - inds[0] + 1));

        inds[0] = 0L;
        inds[1] = nread - 1;
        inds[2] = 1L;

        if (PD_write_alt(out, varname, "integer", idata, 4, inds) == FALSE) {
          fprintf(stderr, "\tWARNING: Could not write '%s'. Ignoring\n", varname);
        }

        delete[] idata;

      } else if ((strcasecmp(typ, "float") == 0) || (strcasecmp(typ, "double") == 0)) {
        // Convert doubles to floats

        float* fdata = new float[varsize];

        // Read the data from the PDB file
        inds[2] = tstride;
        int nread;
        if ((nread = PD_read_as_alt(in, varname, "float", fdata, inds)) == 0) {
          fprintf(stderr, "\tWARNING: Could not read '%s'. Ignoring\n", varname);
          delete[] fdata;
          continue;
        }

        if (nread % (varsize / (inds[1] - inds[0] + 1)) != 0) {
          fprintf(stderr, "ERROR: Accounting error: (%ld, %ld), %d, %d\n", inds[0],
                  inds[1], varsize, nread);
          delete[] fdata;
          continue;
        }

        nread = nread / (varsize / (inds[1] - inds[0] + 1));

        inds[0] = 0L;
        inds[1] = nread - 1;
        inds[2] = 1L;

        if (PD_write_alt(out, varname, "float", fdata, 4, inds) == FALSE) {
          fprintf(stderr, "\tWARNING: Could not write '%s'. Ignoring\n", varname);
        }

        delete[] fdata;

      } else {
        fprintf(stderr, "\tWARNING: '%s' has unrecognised type '%s'. Ignoring\n", varname,
                typ);
      }
    } else {
      // Just copy the data across

      if (strcasecmp(typ, "integer") == 0) {

        int* idata = new int[varsize];

        // Read the data from the PDB file
        if (PD_read_as(in, varname, "integer", idata) == 0) {
          fprintf(stderr, "\t\tWARNING: Could not read variable. Ignoring\n");
          delete[] idata;
          continue;
        }

        if (PD_write_alt(out, varname, "integer", idata, nd, inds) == FALSE) {
          fprintf(stderr, "\tWARNING: Could not write variable '%s'\n", varname);
        }

        delete[] idata;

      } else if ((strcasecmp(typ, "float") == 0) || (strcasecmp(typ, "double") == 0)) {
        // Convert doubles to floats

        float* fdata = new float[varsize];

        // Read the data from the PDB file
        if (PD_read_as(in, varname, "float", fdata) == 0) {
          fprintf(stderr, "\tWARNING: Could not read variable '%s'. Ignoring\n", varname);
          delete[] fdata;
          continue;
        }

        if (PD_write_alt(out, varname, "float", fdata, nd, inds) == FALSE) {
          fprintf(stderr, "\tWARNING: Could not write variable '%s'\n", varname);
        }

        delete[] fdata;
      } else {
        fprintf(stderr, "WARNING: '%s' has unrecognised type '%s'. Ignoring\n", varname,
                typ);
      }
    }
  }

  PD_close(in);
  PD_close(out);

  return 0;
}
