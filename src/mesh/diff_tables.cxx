/*******************************************************************************
 * Lookup tables of functions. Map between names, codes and functions
 *******************************************************************************/


/// Translate between DIFF_METHOD codes, and functions
struct DiffLookup {
  DIFF_METHOD method;
  Mesh::deriv_func func;     // Single-argument differencing function
  Mesh::inner_boundary_deriv_func inner_boundary_func; // Differencing function using forward derivatives
  Mesh::outer_boundary_deriv_func outer_boundary_func; // Differencing function using backward derivatives
  Mesh::upwind_func up_func; // Upwinding function
  Mesh::inner_boundary_upwind_func inner_boundary_up_func; // Upwinding function using forward derivatives
  Mesh::outer_boundary_upwind_func outer_boundary_up_func; // Upwinding function using backward derivatives
};

/// Translate between short names, long names and DIFF_METHOD codes
struct DiffNameLookup {
  DIFF_METHOD method;
  const char* label; // Short name
  const char* name;  // Long name
};

/// Differential function name/code lookup
static DiffNameLookup DiffNameTable[] = { {DIFF_U1, "U1", "First order upwinding"},
					  {DIFF_U2, "U2", "Second order upwinding"},
					  {DIFF_C2, "C2", "Second order central"},
					  {DIFF_W2, "W2", "Second order WENO"},
					  {DIFF_W3, "W3", "Third order WENO"},
					  {DIFF_C4, "C4", "Fourth order central"},
					  {DIFF_U4, "U4", "Fourth order upwinding"},
                      {DIFF_S2, "S2", "Smoothing 2nd order"},
					  {DIFF_FFT, "FFT", "FFT"},
                      {DIFF_NND, "NND", "NND"},
                      {DIFF_SPLIT, "SPLIT", "Split into upwind and central"},
					  {DIFF_DEFAULT}}; // Use to terminate the list

/// First derivative lookup table
static DiffLookup FirstDerivTable[] = { {DIFF_C2, DDX_C2,     DDX_F2, DDX_B2, NULL, NULL, NULL},
					{DIFF_W2, DDX_CWENO2, DDX_F2, DDX_B2, NULL, NULL, NULL},
					{DIFF_W3, DDX_CWENO3, DDX_F4, DDX_B4, NULL, NULL, NULL},
					{DIFF_C4, DDX_C4,     DDX_F4, DDX_B4, NULL, NULL, NULL},
                                        {DIFF_S2, DDX_S2,     NULL,   NULL,   NULL, NULL, NULL},
					{DIFF_FFT, NULL,      NULL,   NULL,   NULL, NULL, NULL},
					{DIFF_DEFAULT}};

/// Second derivative lookup table
static DiffLookup SecondDerivTable[] = { {DIFF_C2, D2DX2_C2, D2DX2_F2, D2DX2_B2, NULL, NULL, NULL},
					 {DIFF_C4, D2DX2_C4, D2DX2_F4, D2DX2_B4, NULL, NULL, NULL},
					 {DIFF_FFT, NULL,    NULL,     NULL,     NULL, NULL, NULL},
					 {DIFF_DEFAULT}};

/// Upwinding functions lookup table
static DiffLookup UpwindTable[] = { {DIFF_U1, NULL, NULL, NULL, VDDX_U1, NULL, NULL},
					{DIFF_U2, NULL, NULL, NULL, VDDX_U2, NULL, NULL}, 
				    {DIFF_C2, NULL, NULL, NULL, VDDX_C2, NULL, NULL},
				    {DIFF_U4, NULL, NULL, NULL, VDDX_U4, NULL, NULL},
				    {DIFF_W3, NULL, NULL, NULL, VDDX_WENO3, NULL, NULL},
				    {DIFF_C4, NULL, NULL, NULL, VDDX_C4, NULL, NULL},
				    {DIFF_DEFAULT}};

/// Flux functions lookup table
static DiffLookup FluxTable[] = { {DIFF_SPLIT, NULL, NULL, NULL, NULL, NULL, NULL},
                                  {DIFF_U1, NULL, NULL, NULL, FDDX_U1, NULL, NULL},
                                  {DIFF_C2, NULL, NULL, NULL, FDDX_C2, NULL, NULL},
                                  {DIFF_C4, NULL, NULL, NULL, FDDX_C4, NULL, NULL},
                                  {DIFF_NND, NULL, NULL, NULL, FDDX_NND, NULL, NULL},
                                  {DIFF_DEFAULT}};

/// First staggered derivative lookup
static DiffLookup FirstStagDerivTable[] = { {DIFF_C2, DDX_C2_stag, DDX_F2_stag, DDX_B2_stag, NULL, NULL, NULL}, 
					    {DIFF_C4, DDX_C4_stag, DDX_F4_stag, DDX_B4_stag, NULL, NULL, NULL},
					    {DIFF_DEFAULT}};

/// Second staggered derivative lookup
static DiffLookup SecondStagDerivTable[] = { {DIFF_C2, D2DX2_C2_stag, D2DX2_F2_stag, D2DX2_B2_stag, NULL, NULL, NULL},
					     {DIFF_DEFAULT}};
                                             //                                             {DIFF_C4, D2DX2_C4_stag, D2DX2_F4_stag, D2DX2_B4_stag, NULL, NULL, NULL},
/// Upwinding staggered lookup
static DiffLookup UpwindStagTable[] = { {DIFF_U1, NULL, NULL, NULL, VDDX_U1_stag, NULL, NULL},
					{DIFF_U2, NULL, NULL, NULL, VDDX_U2_stag, NULL, NULL},
					{DIFF_C2, NULL, NULL, NULL, VDDX_C2_stag, NULL, NULL},
					{DIFF_C4, NULL, NULL, NULL, VDDX_C4_stag, NULL, NULL},
					{DIFF_DEFAULT} };

/// Flux staggered lookup
static DiffLookup FluxStagTable[] = { {DIFF_SPLIT, NULL, NULL, NULL, NULL, NULL, NULL},
                                      {DIFF_U1, NULL, NULL, NULL, FDDX_U1_stag, NULL, NULL},
                                      {DIFF_DEFAULT}};
