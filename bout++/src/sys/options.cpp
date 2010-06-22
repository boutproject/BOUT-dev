/**************************************************************************
 * Reads in the configuration file, supplying
 * an interface to get options
 * 
 * File is an ini file with sections
 * [section]
 * and variables as
 * name = string ; comment
 * 
 * To Do / Known issues
 * ====================
 *
 *  * Currently uses its own hash function. Probably better to just
 *    use an STL map for this class
 *
 * ChangeLog
 * =========
 * 
 * 2010-02-10 Ben Dudson <bd512@york.ac.uk>
 *    
 *    * Adding set methods to allow other means to control code
 *      Intended to help with integration into FACETS
 *
 * 2007-09-01 Ben Dudson <bd512@york.ac.uk>
 *
 *    * Initial version
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
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

#include "globals.h"
#include "options.h"
#include "utils.h"

#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <stdlib.h>

const char SEC_CHAR='_'; ///< Used to separate sections and keys

OptionFile::OptionFile()
{
  noptions = 0;
  def_section = NULL;
}

OptionFile::OptionFile(const char *filename)
{
  noptions = 0;
  read(filename);
  def_section = NULL;
}

OptionFile::~OptionFile()
{
  if(noptions > 0)
    free(option);

  if(def_section != NULL)
    free(def_section);
}


/**************************************************************************
 * Read input file
 **************************************************************************/

int OptionFile::read(const char *format, ...)
{
  FILE *fp;
  int linenr;
  int i, n, p;
  char buffer[MAX_LINELEN];
  char *section; // Current section
  char *s;

  va_list ap;  // List of arguments
  char filename[512];

  if(format == (const char*) NULL) {
    output.write("WARNING: OptionFile read passed NULL filename\n");
    return 1;
  }

  va_start(ap, format);
    vsprintf(filename, format, ap);
  va_end(ap);

  output.write("Reading options file %s\n", filename);

  section = (char*) NULL;

  if((fp = fopen(filename, "rt")) == (FILE*) NULL) {
    output.write("\tOptions file '%s' not found\n", filename);
    return(1);
  }
  
  if((linenr = get_nextline(fp, buffer, MAX_LINELEN, 1)) == -1) {
    output.write("\tEmpty option file '%s'\n", filename);
    return(0);
  }

  do {
    // printf("Line %d:%s:\n", linenr, buffer);
    n = strlen(buffer);
    
    // Check for section
    if(buffer[0] == '[') {
      if(buffer[n-1] != ']') {
	output.write("\t'%s' line %d: Missing ']'\n", filename, linenr);
	return(2);
      }
      buffer[n-1] = 0;
      s = buffer+1;
      
      // Strip lead/trail spaces from the name
      if((n = strip_space(s)) == 0) {
	output.write("\t'%s' line %d: Empty section name\n", filename, linenr);
	return(2);
      }

      // Set section
      section = (char*) realloc(section, n+1);
      strcpy(section, s);
    }else {
      // A key/value pair, separated by a '='
      
      // Find the position of the equality
      p = -1;
      
      for(i=0;i<n;i++) {
	if(buffer[i] == '=') {
	  if(p != -1) {
	    output.write("\t'%s' line %d: Multiple equalities\n", filename, linenr);
	    return(2);
	  }
	  p = i;
	}
      }

      if(p == -1) {
	output.write("\t'%s' line %d: Don't know what this is. Ignoring\n",
		 filename, linenr);
      }else {
      
	s = &(buffer[p+1]); // Now contains value
	buffer[p] = 0;   // Now contains key
      
	if((n = strip_space(buffer)) == 0) {
	  output.write("\t'%s' line %d: Empty key\n", filename, linenr);
	  return(2);
	}
      
	n = strip_space(s);
	
	if( ((s[0] == '\'') || (s[0] == '"')) && (s[n-1] == s[0]) ) {
	  // String in quotes - remove quotes
	  s[n-1] = 0;
	  s++;
	}
      
	add(section, buffer, s, linenr);
      }
    }
  }while((linenr = get_nextline(fp, buffer, MAX_LINELEN, 0)) != -1);

  fclose(fp);
  return(0);
}

int OptionFile::command_line(int argc, char** argv)
{
  output << "Checking command-line options\n";
  // Go through command-line arguments
  for(int i=1;i<argc;i++) {
    // Should contain a "key=value" string
    
    // Find the position of the equality
    int p = -1;
    
    int n = strlen(argv[i]);
    for(int j=0;j<n;j++) {
      if(argv[i][j] == '=') {
	if(p != -1) {
	  // Multiple inequalities - ignore
	  output.write("\tWARNING: Ignoring option %s. Multiple equalities\n");
	  continue;
	}
	p = j;
      }
    }
    
    if(p == -1) {
      /// No equality sign. Assume it's a setting which is set to true
      /// This is so that adding "restart" to the command-line works
      add(NULL, argv[i], "true", -1);
      continue;
    }
    argv[i][p] = 0;
    
    // Add to the hash table
    add(NULL, argv[i], &(argv[i][p+1]), -1);
  }
  return 0;
}

/**************************************************************************
 * Functions to request options
 **************************************************************************/

int OptionFile::getInt(const char *name, int &val)
{
  int i;

  if(name == NULL)
    return 1;

  if((i = find(name)) == -1)
    return(1);

  if(sscanf(option[i].string, "%d", &val) != 1) {
    output.write("\tOption '%s': Integer expected\n", name);
    return(1);
  }

  return(0);
}

int OptionFile::getInt(const char *section, const char *name, int &val)
{
  static char *str;
  static int len = 0;
  int n;
  
  if(name == NULL)
    return 1;

  if(section == NULL)
    return getInt(name, val);

  n = strlen(section) + strlen(name);

  if(n > len) {
    if(len != 0)
      free(str);
    str = (char*) malloc(sizeof(char)*(n+2));
    len = n;
  }
  sprintf(str, "%s%c%s", section, SEC_CHAR, name);
  return( getInt(str, val) );
}

int OptionFile::getReal(const char *name, real &val)
{
  int i;
  double v;
  
  if((i = find(name)) == -1)
    return(1);
  
  if(sscanf(option[i].string, "%lf", &v) != 1) {
    output.write("\tOption '%s': Float expected\n", name);
    return(1);
  }
  
  val = (real) v;

  return(0);
}

int OptionFile::getReal(const char *section, const char *name, real &val)
{
  static char *str;
  static int len = 0;
  int n;
  
  if(name == NULL)
    return 1;

  if(section == NULL)
    return getReal(name, val);

  n = strlen(section) + strlen(name);

  if(n > len) {
    if(len != 0)
      free(str);
    str = (char*) malloc(sizeof(char)*(n+2));
    len = n;
  }
  sprintf(str, "%s%c%s", section, SEC_CHAR, name);
  return( getReal(str, val) );
}

char* OptionFile::getString(const char *name)
{
  if(name == NULL)
    return NULL;

  if(def_section == NULL) {
    int i;
    if((i = find(name)) == -1)
      return((char*) NULL);
    
    return(option[i].string);
  }
  
  return getString(def_section, name);
}

char* OptionFile::getString(const char *section, const char *name)
{
  static char *str;
  static int len = 0;
  int n;

  if(name == NULL)
    return NULL;

  if(section == NULL)
    return getString(name);
  
  n = strlen(section) + strlen(name);

  if(n > len) {
    if(len != 0)
      free(str);
    str = (char*) malloc(sizeof(char)*(n+2));
    len = n;
  }
  sprintf(str, "%s%c%s", section, SEC_CHAR, name);

  int i;
  if((i = find(str)) == -1)
    return((char*) NULL);
  
  return(option[i].string);
}

int OptionFile::getBool(const char *name, bool &val)
{
  int i;
  char c;

  if((i = find(name)) == -1) {
    return(1);
  }

  c = toupper(option[i].string[0]);

  if((c == 'Y') || (c == 'T') || (c == '1')) {
    val = true;
  }else if((c == 'N') || (c == 'F') || (c == '0')) {
    val = false;
  }else {
    output.write("\tOption '%s': Boolean expected\n", name);
    return(1);
  }

  return(0);
}

int OptionFile::getBool(const char *section, const char *name, bool &val)
{
  static char *str;
  static int len = 0;
  int n;
  
  if(name == NULL)
    return 1;

  if(section == NULL)
    return getBool(name, val);

  n = strlen(section) + strlen(name);

  if(n > len) {
    if(len != 0)
      free(str);
    str = (char*) malloc(sizeof(char)*(n+2));
    len = n;
  }
  sprintf(str, "%s%c%s", section, SEC_CHAR, name);
  return( getBool(str, val) );
}

/**************************************************************************
 * Set methods
 *
 * Convert given options to strings
 **************************************************************************/

int OptionFile::set(const char *name, int val)
{
  char buffer[128];
  
  snprintf(buffer, 127, "%d", val);
  return set(name, buffer);
}

int OptionFile::set(const char *name, real val)
{
  char buffer[128];
  
  snprintf(buffer, 127, "%le", val);
  return set(name, buffer);
}

int OptionFile::set(const char *name, bool val)
{
  if(val) {
    return set(name, "true");
  }
  
  return set(name, "false");
}

int OptionFile::set(const char *name, const char *str)
{
  int n = strlen(name)+1;
  if(def_section != (char*) NULL)
    n += strlen(def_section) + 1;

  char* s = (char*) malloc(n);
  if(def_section != (char*) NULL) {
    sprintf(s, "%s%c%s", def_section, SEC_CHAR, name);
  }else
    strcpy(s, name);

  // See if the name is already set
  n = find(s);
  if(n != -1) {
    output.write("WARNING: set method changing option '%s' from '%s' to '%s'\n", 
		 s, option[n].string, str);
    
    free(option[n].string);
    free(s);
  }else {
    // Allocate a new option
    if(noptions == 0) {
      option = (t_option*) malloc(sizeof(t_option));
    }else {
      option = (t_option*) realloc(option, (noptions+1)*sizeof(t_option));
    }
    n = noptions;
    noptions++;
    
    option[n].name = s;
    option[n].hash = hash_string(name);
  }
  
  option[n].string = copy_string(str);
  
  return 0;
}

/**************************************************************************
 * Private functions
 **************************************************************************/

void OptionFile::add(const char *section, const char *name, const char *string, int linenr)
{
  int n;

  char *s;

  if(name == NULL)
    return;

  n = strlen(name)+1;
  if(section != (char*) NULL)
    n += strlen(section) + 1;
  s = (char*) malloc(n);
  if(section != (char*) NULL) {
    sprintf(s, "%s%c%s", section, SEC_CHAR, name);
  }else
    strcpy(s, name);

  int id = find(s);
  if(id != -1) {
    if(linenr == -1) {
      // On the command line - override
      output.write("\tWARNING: Option '%s' over-ridden on command-line to '%s'\n", s, string);
      // Free the old strings
      free(option[id].name);
      free(option[id].string);
    }else {
      output.write("\tWARNING: Variable '%s' re-defined on line %d. Keeping first occurrence\n", s, linenr);
      free(s);
      return;
    }
  }else {
    // Allocate memory
    
    if(noptions == 0) {
      option = (t_option*) malloc(sizeof(t_option));
    }else {
      option = (t_option*) realloc(option, (noptions+1)*sizeof(t_option));
    }
    
    id = noptions;
    noptions++;
  }
    
  // Copy name across

  option[id].name = s;

  // Calculate hash

  option[id].hash = hash_string(s);

  // Copy string across
  option[id].string = copy_string(string);
}

int OptionFile::find(const char *name)
{
  int i, n;
  int hash;

  if(name == NULL)
    return -1;

  hash = hash_string(name);

  n = -1;

  for(i=0;i<noptions;i++) {
    if(option[i].hash == hash) {
      // Compare strings to be sure
      if(strcasecmp(option[i].name, name) == 0) {
	n = i;
	break;
      }
    }
  }
  
  return(n);
}

unsigned int OptionFile::hash_string(const char *string)
{
  int i, n;
  unsigned int hash;

  if(string == NULL)
    return 0;

  n = strlen(string);
  hash = 0;

  for(i=0;i<n;i++) {
    hash += (unsigned int) toupper(string[i]);
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;

  return(hash);
}

// Strips leading and trailing spaces from a string
int OptionFile::strip_space(char *string)
{
  int s, i, n;
  
  if(string == NULL)
    return 0;

  n = strlen(string);

  if(n == 0)
    return(0);

  // Strip leading space
  
  s = 0;
  while(isspace(string[s]) && (string[s] != 0))
    s++;
  
  if(string[s] == 0) {
    string[0] = 0;
    return(0);
  }

  // Strip trailing space

  i = n-1;
  while(isspace(string[i]))
    i--;

  n = i+1 - s; // Length of new string
  for(i=0;i<n;i++)
    string[i] = string[i+s];
  string[n] = 0;

  return(n);
}

int OptionFile::get_nextline(FILE *fp, char *buffer, int maxbuffer, int first)
{
  // Returns the next useful line, stripped of comments and whitespace

  int i, n;
  static int linenr = 0;
  
  if(first)
    linenr = 0;

  do{
    buffer[0] = 0;
    /* Get a line from the file */
    fgets(buffer, maxbuffer-1, fp);
    buffer[maxbuffer-1] = 0; // ensure always have terminating zero
    linenr++;

    n = strlen(buffer);

    // Strip out comments
    for(i=0;i<n;i++) {
      if((buffer[i] == COMMENT_CHAR) || (buffer[i] == '#'))
	buffer[i] = 0;
    }
    
    // Strip out whitespace
    n = strip_space(buffer);

  }while((feof(fp) == 0) && (n == 0));

  if(n == 0)
    linenr = -1;

  return(linenr);
}


// New interface

void OptionFile::setSection(const char *name) // Set the default section
{
  int n;

  if(def_section != NULL)
    free(def_section);

  if(name == NULL) {
    def_section = NULL;
    return;
  }

  n = strlen(name);
  if(n == 0) {
    def_section = NULL;
    return;
  }

  def_section = (char*) malloc(n+1);
  
  strcpy(def_section, name);
}
 
int OptionFile::get(const char *name, int &val, const int def)
{
  if(name == NULL) {
    output.write("WARNING: NULL integer option requested\n");
    return 1;
  }
    
  if(def_section == NULL) {
    if(getInt(name, val)) {
      val = def;
      output.write("\tOption %s = %d (default)\n", name, def);
      return 1;
    }else {
      output.write("\tOption %s = %d\n", name, val);
    }
  }else {
    if(getInt(def_section, name, val)) {
      val = def;
      output.write("\tOption %s / %s = %d (default)\n", def_section, name, def);
      return 1;
    }else {
      output.write("\tOption %s / %s = %d\n", def_section, name, val);
    }
  }
  return 0;
}

int OptionFile::get(const char *name, real &val, const real def)
{
  if(name == NULL) {
    output.write("WARNING: NULL real option requested\n");
    return 1;
  }
  
  if(def_section == NULL) {
    if(getReal(name, val)) {
      val = def;
      output.write("\tOption %s = %e (default)\n", name, def);
      return 1;
    }else {
      output.write("\tOption %s = %e\n", name, val);
    }
  }else {
    if(getReal(def_section, name, val)) {
      val = def;
      output.write("\tOption %s / %s = %e (default)\n", def_section, name, def);
      return 1;
    }else {
      output.write("\tOption %s / %s = %e\n", def_section, name, val);
    }
  }
  return 0;
}

int OptionFile::get(const char *name, bool &val, const bool def)
{
  if(name == NULL) {
    output.write("WARNING: NULL boolean option requested\n");
    return 1;
  }

  if(def_section == NULL) {
    if(getBool(name, val)) {
      val = def;
      if(def) {
	output.write("\tOption %s = true (default)\n", name);
      }else
	output.write("\tOption %s = false (default)\n", name);
      return 1;
    }else {
      if(val) {
	output.write("\tOption %s = true\n", name);
      }else
	output.write("\tOption %s = false\n", name);
    }
  }else {
    if(getBool(def_section, name, val)) {
      val = def;
      if(def) {
	output.write("\tOption %s / %s = true (default)\n", def_section, name);
      }else
	output.write("\tOption %s / %s = false (default)\n", def_section, name);
      return 1;
    }else {
      if(val) {
	output.write("\tOption %s / %s = true\n", def_section, name);
      }else
	output.write("\tOption %s / %s = false\n", def_section, name);
    }
  }
  return 0;
}

int OptionFile::get(const char *section, const char *name, int &val, const int def)
{
  if(name == NULL) {
    output.write("WARNING: NULL integer option requested\n");
    return 1;
  }
  
  if(section == NULL)
    return get(name, val, def);
  
  if(getInt(section, name, val)) {
    val = def;
    output.write("\tOption %s / %s = %d (default)\n", section, name, def);
    return 1;
  }else {
    output.write("\tOption %s / %s = %d\n", section, name, val);
  }
  return 0;
}

int OptionFile::get(const char *section, const char *name, real &val, const real def)
{
  if(name == NULL) {
    output.write("WARNING: NULL real option requested\n");
    return 1;
  }
  
  if(section == NULL)
    return get(name, val, def);

  if(getReal(section, name, val)) {
    val = def;
    output.write("\tOption %s / %s = %e (default)\n", section, name, def);
    return 1;
  }else {
    output.write("\tOption %s / %s = %e\n", section, name, val);
  }
  return 0;
}

int OptionFile::get(const char *section, const char *name, bool &val, const bool def)
{
  if(name == NULL) {
    output.write("WARNING: NULL boolean option requested\n");
    return 1;
  }
  
  if(section == NULL)
    return get(name, val, def);
  
  if(getBool(section, name, val)) {
    val = def;
    if(def) {
      output.write("\tOption %s / %s = true (default)\n", section, name);
    }else
      output.write("\tOption %s / %s = false (default)\n", section, name);
    return 1;
  }else {
    if(val) {
      output.write("\tOption %s / %s = true\n", section, name);
    }else
      output.write("\tOption %s / %s = false\n", section, name);
  }
  return 0;
}

