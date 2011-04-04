; Compile shared library using IDL's routines
; 
; Set the environment variable PACT to point to the PACT library

PACT=GETENV("PACT")

MAKE_DLL, ["pdb2idl", "utils"], "pdb2idl", $
  ["pdb_open", "pdb_close", "idl_legalvar", $
   "pdb_list_vars", "pdb_list_varname", "pdb_list_free", $
   "pdb_query_type", $
   "pdb_query_ndim", "pdb_query_dims", $
   "pdb_read_string", $
   "pdb_read_0d_int", "pdb_read_1d_int", $
   "pdb_read_0d_float", "pdb_read_1d_float", "pdb_read_2d_float", "pdb_read_3d_float", "pdb_read_4d_float", $
   "pdb_write_string", $
   "pdb_write_0d_int", "pdb_write_1d_int", $
   "pdb_write_0d_float", "pdb_write_1d_float", "pdb_write_2d_float", "pdb_write_3d_float", "pdb_write_4d_float", $
   "pdb_write_0d_double", "pdb_write_1d_double"], $
  COMPILE_DIR=".", $
  EXTRA_CFLAGS="-I"+PACT+"/include", $
  EXTRA_LFLAGS="-L"+PACT+"/lib -lpdb -lpml -lscore", $
  DLL_PATH=DLL_PATH

PRINT, "Shared library is ", DLL_PATH

exit
