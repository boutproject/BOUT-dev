from common import *

# read tables
func_tables={}
first_entry={}
cfunc=None
with open("tables_cleaned.cxx","r") as f:
    inFunc=0;
    for line in f:
        inFunc+=line.count('{')
        if inFunc:
            if cfunc is None:
                cfunc=[]
            cfunc.append(line);
        inFunc-=line.count('}')
        if not inFunc:
            if not cfunc is None:
                name=cfunc[0].split(" ")[2].split("[")[0]
                if len(name) > 2:
                    #print name
                    func_tables[name]={}
                    for e in cfunc:
                        e=e.split('{')[-1]
                        e=e.split('}')[0]
                        e=e.split(',')
                        en=[]
                        for diff in e:
                            en.append(diff.strip())
                        cn=en[0]
                        if cn=='DIFF_W3':
                            continue
                        #print cn
                        if not name in first_entry:
                            first_entry[name]=cn
                        if en[1:7] == ['NULL']*6:
                            continue
                        func_tables[name][cn]=en[1:7]
                cfunc=None
descriptions=func_tables.pop("DiffNameTable")

# fields=['Field3D','Field2D']
# dirs={}
# dirs['Field3D']=['x','y','z']
# dirs['Field2D']=['x','y']
funcname={ 'FirstDerivTable' : 'indexDD%s',
           'SecondDerivTable' : 'indexD2D%s2',
           'FirstStagDerivTable' : 'indexDD%s',
           'SecondStagDerivTable' : 'indexD2D%s2'}

funcs_to_gen=[]
for t in func_tables:
    #print t
    #if t== "FirstDerivTable":
    #    print func_tables[t]
    func_tables[t].pop('DIFF_DEFAULT')
    fu=func_tables[t].itervalues().next()
    if fu[1] != "NULL": # not a flux/upwind scheeme
        #print t
        if fu[1][-4:]== "stag":
            stag=True
        else:
            stag=False
        #print t,"is %s staggered"%stag
        for field in fields:
            for d in dirs[field]:
                warn()
                if not stag:
                    myname=funcname[t]%d.upper()+"_non_stag"
                else:
                    myname=funcname[t]%d.upper()+"_stag"
                print "const",field, myname , "(const",field,"&f, CELL_LOC outloc, DIFF_METHOD method) {"
                print "  if (method == DIFF_DEFAULT){"
                print "    method = default_%s;"%t[:-5]# drop 'Table'
                print "  }"
                print "  if (outloc == CELL_DEFAULT){"
                print "    outloc = f.getLocation();"
                print "  }"
                print "  switch (method){"
                for method in func_tables[t]:
                    print "  case",method,":"
                    if stag:
                        print "    if (outloc == CELL_%sLOW){"%d.upper()
                        print "      if (f.getLocation() == CELL_CENTRE){"
                        print "        return %s_on_%s(f);"%(funcname[t]%d.upper(),method) #func_tables[t][method][0]
                        print "      } else {"
                        print "        return %s_on_%s(interp_to(f,CELL_CENTRE));"%(funcname[t]%d.upper(),method)
                        print "      };"
                        print "    } else {" # inloc must be CELL_%sLOW
                        print "      if (outloc == CELL_CENTRE){"
                        print "        return %s_off_%s(f);"%(funcname[t]%d.upper(),method)
                        print "      } else {"
                        print "        return interp_to(%s_off_%s(f),outloc);"%(funcname[t]%d.upper(),method)
                        print "      };"
                        print "    }"
                        stags=['on','off']
                    else:
                        print "    if (outloc == f.getLocation()){"
                        print "      return %s_norm_%s(f);"%(funcname[t]%d.upper(),method)
                        print "    } else {"
                        print "      return interp_to(%s_norm_%s(f),outloc);"%(funcname[t]%d.upper(),method)
                        print "    }"
                        stags=['norm']
                    for mstag in stags:
                        funcs=func_tables[t][method]
                        forward=funcs[1]
                        if forward=='NULL':
                            forward=func_tables[t][first_entry[t]][1]
                        backward=funcs[2]
                        if backward=='NULL':
                            backward=func_tables[t][first_entry[t]][2]
                        funcs_to_gen.append(["%s_%s_%s"%(funcname[t]%d.upper(),mstag,method),field,d,mstag,funcs[0],forward,backward])
                    print "    break;"
                    #print "    }"
                print "  default:"
                print "    throw BoutException(\"%s CartesianMesh::"%field +myname,'unknown method %d.\\nNote FFTs are not (yet) supported.",method);'
                print "  }; // end switch"
                print "};"
                print
                
    #else:
    #    print fu
    #for meth in func_tables[t]:
    #    print meth
    #print func_tables[t]
#print first_entry

headers=""
for func in ["indexDD%s", "indexD2D%s2"]:
    for field in fields:
        for d in dirs[field]:
            warn()
            function_header="  virtual const "+field+" "+func%d.upper()
            function_header+="(const "+field+" &f, CELL_LOC outloc, DIFF_METHOD method)";
            headers+=function_header+";\n"
            function_header="const "+field+" CartesianMesh::"+func%d.upper()
            function_header+="(const "+field+" &f, CELL_LOC outloc, DIFF_METHOD method)";
            print function_header," {"
            print "  if ((outloc == CELL_%sLOW) != (f.getLocation() == CELL_%sLOW)){"%(d.upper(),d.upper())
            print "    // we are going onto a staggered grid or coming from one"
            print "    return",func%d.upper()+"_stag","(f,outloc,method);"
            print "  } else {"
            print "    return",func%d.upper()+"_non_stag","(f,outloc,method);"
            print "  }"
            print "}"
            print

with open("generated_header.hxx","w") as f:
    f.write(headers)


import sys
#funcs_to_gen=funcs_to_gen[0:2]
sys.stdout=open("generated_stencils.cxx","w")
from gen_stencils import gen_functions_normal
gen_functions_normal(funcs_to_gen)
