from common import *
import sys

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
                        #NI not implemented :
                        if cn=='DIFF_W3':
                            continue
                        if cn=='DIFF_SPLIT':
                            continue
                        #print cn
                        if not name in first_entry:
                            first_entry[name]=cn
                        if en[1:8] == ['NULL']*7:
                            continue
                        func_tables[name][cn]=en[1:8]
                cfunc=None
descriptions=func_tables.pop("DiffNameTable")
descriptions.pop("DIFF_DEFAULT")

# import sys
# print >> sys.stderr , descriptions
# exit(1)

funcname={ 'FirstDerivTable'      : 'indexDD%s',
           'FirstStagDerivTable'  : 'indexDD%s',
           'SecondDerivTable'     : 'indexD2D%s2',
           'SecondStagDerivTable' : 'indexD2D%s2',
           'UpwindTable'          : 'indexVDD%s',
           'UpwindStagTable'      : 'indexVDD%s',
           'FluxTable'            : 'indexFDD%s',
           'FluxStagTable'        : 'indexFDD%s'
}

funcs_to_gen=[]
default_methods=dict()
duplicates(func_tables.keys())
#import sys
#for f in func_tables:
#    fu=func_tables[f]
#    for k in fu:
#        print >>sys.stderr, k,fu[k]
for t in func_tables:
    try:
        func_tables[t].pop('DIFF_DEFAULT')
    except:
        #print >>sys.stderr,t
        #print >>sys.stderr,func_tables[t]
        pass
    try:
        func_tables[t].pop('DIFF_SPLIT')
    except:
        pass
    fu=func_tables[t].itervalues().next()
    if fu[1] != "NULL": # not a flux/upwind scheeme
        flux=False
        upwind=False
    else:
        if fu[3] != "NULL":
            upwind=True
        else:
            upwind=False
        flux=True
    if upwind:
        pos=3
    elif flux:
        pos=4
    else:
        pos=0
    if fu[pos][-4:]== "stag":
        stag=True
    else:
        stag=False
    if True:
        #for t in func_tables:
        #print >>sys.stderr , t #func_tables
        duplicates(fields)
        for field in fields:
            duplicates(dirs[field])
            for d in dirs[field]:
                warn()
                try:
                    if not stag:
                        myname=funcname[t]%d.upper()+"_non_stag"
                    else:
                        myname=funcname[t]%d.upper()+"_stag"
                except:
                    import sys
                    print >>sys.stderr,funcname[t]
                    exit(3)
                if flux:
                    inp="(const "+field+" &v, const "+field+" &f, "
                else:
                    inp="(const "+field+" &f, "
                print "const",field, myname , inp,"CELL_LOC outloc, DIFF_METHOD method) {"
                print "  if (method == DIFF_DEFAULT){"
                print "    method = default_%s_%s;"%(d,t[:-5]+("Deriv" if flux else ""))# drop 'Table'
                print "  }"
                print "  if (outloc == CELL_DEFAULT){"
                print "    outloc = f.getLocation();"
                print "  }"
                print "  switch (method){"
                default_methods["default_%s_%s"%(d,t[:-5])]=func_tables[t]
                duplicates(func_tables[t].keys())
                for method in func_tables[t]:
                    print "  case",method,":"
                    if flux:
                        f="v,f"
                    else:
                        f="f"
                    if flux:
                        # f.getLocation() == outloc guaranteed
                        if stag:
                            print "    if (outloc == CELL_%sLOW){"%d.upper()
                            print "      return %s_on_%s(interp_to(v,CELL_CENTRE),f);"% \
                                (funcname[t]%d.upper(),method)
                            print "    } else {" # inloc must be CELL_%sLOW
                            print "      return interp_to(%s_off_%s(v,interp_to(f,CELL_CENTRE)),outloc);"% \
                                (funcname[t]%d.upper(),method)
                            print "    }"
                            stags=['on','off']
                        else:
                            print "    if (v.getLocation() == f.getLocation()){"
                            print "      return interp_to(%s_norm_%s(v,f),outloc);"%(funcname[t]%d.upper(),method)
                            print "    } else {"
                            print "      return interp_to(%s_norm_%s(interp_to(v,CELL_CENTRE),interp_to(f,CELL_CENTRE)),outloc);"%(funcname[t]%d.upper(),method)
                            print "    }"
                            stags=['norm']
                    else: # not flux
                        if stag:
                            print "    if (outloc == CELL_%sLOW){"%d.upper()
                            print "      return %s_on_%s(interp_to(%s,CELL_CENTRE));"%(funcname[t]%d.upper(),method,f)
                            print "    } else {" # inloc must be CELL_%sLOW
                            print "      return interp_to(%s_off_%s(%s),outloc);"% \
                                (funcname[t]%d.upper(),method,f)
                            print "    }"
                            stags=['on','off']
                        else:
                            print "    return interp_to(%s_norm_%s(%s),outloc);"%(funcname[t]%d.upper(),method,f)
                            stags=['norm']
                    for mstag in stags:
                        funcs=func_tables[t][method]
                        if funcs[0]=='NULL':
                            #print >> sys.stderr, funcs
                            funcs[1:3]=funcs[5:7]
                            funcs[0]=funcs[3]
                            if funcs[0]=='NULL':
                                funcs[0]=funcs[4]
                            #print >> sys.stderr, funcs
                        forward=funcs[1]
                        if forward=='NULL':
                            try:
                                forward=func_tables[t][first_entry[t]][1]
                                #if forward=='NULL':
                                #    raise
                            except:
                                print >> sys.stderr,t
                                print >> sys.stderr,first_entry
                                print >> sys.stderr,func_tables[t]
                                raise
                        backward=funcs[2]
                        if backward=='NULL':
                            backward=func_tables[t][first_entry[t]][2]
                        #print >> sys.stderr, funcname[t] , d, mstag, method, field
                        funcs_to_gen.append(["%s_%s_%s"%(funcname[t]%d.upper(),mstag,method),field,d,mstag,funcs[0],forward,backward,flux])
                    print "    break;"
                    #print "    }"
                print "  default:"
                print "    throw BoutException(\"%s AiolosMesh::"%field +myname,'unknown method %d.\\nNote FFTs are not (yet) supported.",method);'
                print "  }; // end switch"
                print "}"
                print

    #else:
    #    print fu
    #for meth in func_tables[t]:
    #    print meth
    #print func_tables[t]
#print first_entry

headers=""
for func in ["indexDD%s", "indexD2D%s2","indexVDD%s","indexFDD%s"]:
    flux=True
    if func.find("indexD") > -1:
        flux=False
    for field in fields:
        for d in dirs[field]:
            warn()
            sig="("
            if flux:
                sig+="const "+field+" &v,"
            sig += "const "+field+" &f"
            if field=="Field3D" or flux:
                sig+=", CELL_LOC outloc, DIFF_METHOD method";
            if func%d.upper() in ["indexDDZ", "indexD2DZ2"]:
                sig+=",bool ignored";
            sig+=")"
            function_header="  virtual const "+field+" "+func%d.upper()
            function_header+=sig
            if  not (field == "Field3D" and func[5]=='V'):
                #print >> sys.stderr , func%d    #     indexD2DZ2
                #if not (field == "Field3D" and func%d.upper() == "indexD2DZ2"):
                function_header+=" override;\n"
                #else:
                #    function_header+=";\n"
            else:
                function_header+=""";
virtual const Field3D indexVDD%s(const Field &v,const Field &f, CELL_LOC outloc, DIFF_METHOD method) override{
  return indexVDD%s(dynamic_cast<const Field3D &>(v),dynamic_cast<const Field3D &>(f),outloc,method);
}"""%(d.upper(),d.upper())
            headers+=function_header
            function_header="const "+field+" AiolosMesh::"+func%d.upper()
            function_header+=sig
            if flux:
                f="v, f"
            else:
                f="f"
            print function_header," {"
            if field != "Field3D" and not flux:
                print "  CELL_LOC outloc=CELL_DEFAULT;"
                print "  DIFF_METHOD method=DIFF_DEFAULT;"
            print "  if (outloc == CELL_DEFAULT) {"
            print "    outloc=f.getLocation();"
            print "  }"
            if flux:
                print "  if (outloc != f.getLocation()) {"
                print '    throw BoutException("AiolosMesh::index?DDX: Unhandled case for shifting.\\n\
f.getLocation()==outloc is required!");'
                print "  }"
            #print '  output.write("Using aiolos mesh for %s\\n");'%(func%d.upper())
            print "  if ((outloc == CELL_%sLOW) != (f.getLocation() == CELL_%sLOW)){"% \
                (d.upper(),d.upper())
            print "    // we are going onto a staggered grid or coming from one"
            print "    return",func%d.upper()+"_stag","("+f+",outloc,method);"
            print "  } else {"
            print "    return",func%d.upper()+"_non_stag","("+f+",outloc,method);"
            print "  }"
            print "}"
            print

with open("generated_header.hxx","w") as f:
    f.write(headers)


import sys
#funcs_to_gen=funcs_to_gen[0:2]
tmp=[]
for fu in funcs_to_gen:
    tmp.append(fu[0]+fu[1])
duplicates(tmp)
#seen = set()
#uniq = []
# for x in funcs_to_gen:
#     xs=x[0]+x[1]
#     if xs not in seen:
#         uniq.append(x)
#         seen.add(xs)
# funcs_to_gen=uniq
guards_=[]
sys.stdout=open("generated_stencils.cxx","w")
from gen_stencils import gen_functions_normal
gen_functions_normal(funcs_to_gen)

sys.stdout=open("generated_init.cxx","w")

descriptions_cleaned=dict()
for d in descriptions:
    #print >> sys.stderr, d
    descriptions_cleaned[d]=descriptions[d][1].strip('"');
for d in dirs['Field3D']:
    warn()
    for i in ['First','Second','Upwind','Flux']:
        stags=['','Stag']
        if i in ['First','Second']:
            table="DerivTable"
        else:
            table="Table"
            #stags=[""]
        for stag in stags:
            print 'DIFF_METHOD default_%s_%s%sDeriv;'%(d,i,stag)

warn()
print "void AiolosMesh::derivs_init(Options * option) {"
print "  std::string name;"
print "  Options * dirOption;"
print '  Options * defOption = option->getSection("diff");'
for d in dirs['Field3D']:
    print "  output.write(\"\\tSetting derivatives for direction %s:\\n\");"%d
    print '  dirOption = option->getSection("dd%s");'%d
    print
    for i in ['First','Second','Upwind','Flux']:
        stags=['','Stag']
        if i in ['First','Second']:
            table="DerivTable"
        else:
            table="Table"
            #stags=[""]
        for stag in stags:
            warn()
            print '  // Setting derivatives for dd%s and %s'%(d,i+stag)
            print ' ',
            for option in ['dirOption','defOption']:
                if stag == "Stag":
                    names = [ i+stag,i,"all" ]
                else:
                    names = [ i,"all" ]
                for name in names:
                    print 'if (%s->isSet("%s")){'%(option,name)
                    print '    %s->get("%s",name,"");'%(option,name)
                    print '  } else',
            print '{'
            if i == 'Second' and stag == "Stag":
                print '    name="C2";'
            elif i == 'Flux':
                print '    name="SPLIT";'
            elif i == 'Upwind':
                print '    name="U1";'
            else:
                print '    name="C2";'
            #elif i == ''
            print '  }'
            print ' ',
            options=""
            for avail in func_tables[i+stag+table]:
                print 'if (strcasecmp(name.c_str(),"%s")==0) {'%avail[5:]
                print '    default_%s_%s%sDeriv = %s;'%(d,i,stag,avail)
                print '    output.write("\t%15s : %s\\n");'%(i+stag,descriptions_cleaned[avail]);
                print '  } else',
                options+="\\n * %s"%avail[5:]
            print '{'
            print '    throw BoutException("Dont\'t know what diff method to use for %s (direction %s, tried to use %s)!\\nOptions are:%s",name.c_str());'%(i+stag,d,'%s',options)
            print '  }'
print "}"
#exit(1)
