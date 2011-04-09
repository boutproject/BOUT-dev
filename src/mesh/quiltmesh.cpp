#include "quiltmesh.h"
#include "bout_exception.h"

#include <algorithm>
#include <numeric>
using std::accumulate;
#include <cmath>

QuiltMesh::~QuiltMesh()
{

}

int QuiltMesh::load()
{
  // Get number of processors
  int NPES, MYPE;
  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);
  
  // Get number of regions
  int nregions;
  if(get(nregions, "nregions"))
    throw BoutException("Mesh file doesn't have nregions variable. Incorrect format?\n");
  
  if(nregions > NPES)
    throw BoutException("Can't divide %d regions between %d processors\n", nregions, NPES);
  
  // Get array of grid sizes
  vector<int> nx = readInts("nx", nregions);
  vector<int> ny = readInts("ny", nregions);
  
  ///////////////////////////////////////////////////////////
  // Partition this grid between NPES processors.
  // Try to make the number of points on each processor approximately
  // equal, and make each domain as square as possible (minimise comms)
  
  // Number of points in each region ntot = nx * ny
  vector<int> ntot(nregions);
  transform(nx.begin(), nx.end(), 
            ny.begin(),
            ntot.begin(),
            std::multiplies<int>());
  
  // Sum ntot to get total number of points in the grid
  int alln =  accumulate( ntot.begin(), ntot.end(), 0 );
  
  // Allocate processors to regions
  vector<int> nproc(nregions);
  for(int i=0;i<nregions;i++) {
    nproc[i] = 1; // Need at least one per region
    // Share the others out
    nproc[i] += (int) ((NPES-nregions) * ntot[i]) / alln;
  }
  // Find number of processors remaining
  int npremain = NPES - accumulate( nproc.begin(), nproc.end(), 0 );
  
  if(npremain != 0) {
    // Assign extra points to the regions with highest load
    vector<BoutReal> load(nregions);
    for(int i=0;i<nregions;i++)
      load[i] = ((BoutReal) ntot[i]) / ((BoutReal) nproc[i]);
    
    for(int p=0;p<npremain;p++) {
      // Get maximum load index
      int ind = 0;
      BoutReal maxload = load[0];
      for(int i=1;i<nregions;i++)
        if(load[i] > maxload) {
          ind = i;
          maxload = load[i];
        }else if(fabs(load[i] - maxload) < 1.e-5) {
          // Equal loads, so add to the one with the largest number of points
          // (minimise difference in loads after)
          if(ntot[i] > ntot[ind])
            ind = i;
        }
      // Add a processor to this region
      nproc[ind]++;
      load[ind] = ((BoutReal) ntot[ind]) / ((BoutReal) nproc[ind]);
    }
  }
  
  ///////////////////////////////////////////////////////////
  // Each region now has ntot points and nproc processors
  // Divide up each region.
  int proc = 0; // Processor number
  for(int r = 0; r < nregion; r++) { // Loop over the regions
    // Calculate number of processors in X for a square domain (mxsub = mysub)
    int nxp_sq = (int) ((BoutReal) nx[r])*sqrt(((BoutReal) nproc[r]) / ((BoutReal) ntot[r]));
    
    for(int nxp=nxp_sq; nxp > 0; nxp--) {
      if(nproc[r] % nxp == 0) {
        break;
      }
    }
    int nyp = nproc[r] / nxp;

    // Divide up into nyp strips
    int nyall = (int) ny[r] / nyp;
    int nyextra = ny[r] % nyp;
    
    int nxall = (int) nx[r] / nxp;
    int nxextra = nx[r] % nxp;

    int x0 = 0;
    for(int xp = 0; xp < nxp; xp++) {
      int nxsub = nxall; // nx on these processors
      if(xp < nxextra)
        nxsub++;
      for(int yp = 0; yp < nyp; yp++) {
        int y0 = 0;
        int nysub = nyall; // ny on these processors
        if(yp < nyextra)
          nysub++;

        // Create a new domain
        MeshDomain *d = new MeshDomain;
        d->region = r;
        d->nx = nxsub;
        d->ny = nysub;
        d->x0 = x0;
        d->y0 = y0;
        d->proc = proc;
        if(proc == MYPE)
          mydomain = d;
        
        
        domain.push_back(d); // Add this domain
        
        proc++;
        y0 += nysub;
      }
    }
    x0 += nxsub;
  }
}

const vector<int> QuiltMesh::readInts(const string &name, int n)
{
  vector<int> result;
  
  // First get a data source
  GridDataSource* s = findSource(name);
  if(s) {
    s->open(name);
    s->setOrigin();
    if(!s->fetch(result, name, n)) {
      // Error reading
      s->close();
      throw BoutException("Could not read integer array '%s'\n", name.c_str());
    }
    s->close();
  }else {
    // Not found
    throw BoutException("Missing integer array %s\n", name.c_str());
  }
  
  return result;
}
