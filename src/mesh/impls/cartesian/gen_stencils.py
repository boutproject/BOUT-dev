# Define some constants, functions
def warn():
    print "// This file auto-generated - do not edit!"

off_diff=dict()
off_diff['on']={
    'ppp':'ppp',
    'pp' :'pp',
    'p'  :'p',
    'm'  :'c',
    'mm' :'m',
    'mmm':'mm'}
off_diff['off']={
    'ppp':'pp',
    'pp' :'p',
    'p'  :'c',
    'm'  :'m',
    'mm' :'mm',
    'mmm':'mmm'}
off_diff['norm']={
    'ppp':'ppp',
    'pp' :'pp',
    'p'  :'p',
    'c'  :'c',
    'm'  :'m',
    'mm' :'mm',
    'mmm':'mmm'}

# Read the stencils
functions=[]
inFunc=0
with open("stencils_cleaned.cxx","r") as f:
    fi=-1
    for line in f:
        line=line.strip()
        line=line.replace("result","result_") # avoid to overwrite result
        if line == '':
            continue
        if line.find("{") > -1: #start of function
            if inFunc==0:
                functions.append([])
                fi+=1
            inFunc+=1
        if inFunc:
            functions[fi].append(line);
        if line.find("}") > -1:
            inFunc-=1

func_db=dict()

dir=['x', 'y', 'z']
#dir=['x','y']
for d in dir:
    func_db[d]=[]

class function:
    name=""
    body=[]
    stag=False
    flux=False
    

functions_=functions
functions=[]
for f in functions_:
    func=function()
    func.body=f
    if f[0].find("BoutReal") == 0 and f[0].find('stencil') > 0:
        #for i in f:
        #    print i
        func.name=f[0].split()[1].split("(")[0]
        if func.name=='DDX_CWENO3':
            #print << sys.stderr , "Skipping DDX_CWENO3 ..."
            continue;
        #print fname
        if (f[0].find("BoutReal V") == -1 and f[0].find("BoutReal F") == -1 ):
            # They have one stencil
            #print f[0]
            func.flux=False
        else:
            func.flux=True
        if f[0].find("stag") > 0:
            # They are staggered
            func.stag=True
        else:
            func.stag=False
        functions.append(func)
for f in functions:
    if f.flux == False:
        if f.stag:
            for d in dir:
                tmp="cart_diff_%s_%s_%%s"%(d,f.name)
                func_db[d].append([f.name,"NULL",tmp%"off",tmp%"on"])
        else:
            #not staggerd
            for d in dir:
                tmp="cart_diff_%s_%s_%%s"%(d,f.name)
                func_db[d].append([f.name,"cart_diff_%s_%s_norm"%(d,f.name),"NULL","NULL"])
for db in func_db:
    func_db[db].append(["NULL","NULL","NULL","NULL"])
                
for f in functions:
    if f.flux == False:
            if f.stag:
                modes=['on','off']
            else:
                modes=['norm']
            for d in dir:
                tmp="cart_diff_%s_%s_%%s"%(d,f.name)
                for mode in modes:
                        warn()
                        print """Field3D %s(const Field3D in){
  Field3D result;
  //result.allocate();
  result=0;"""%(tmp%mode)
                        if d != 'z':
                            print """  for (auto i: result.region(RGN_NO%s)){"""%d.upper()
                        else:
                            print """  for (auto i: result){"""
                        for line in f.body[1:]:
                            ret="return"
                            pos=line.find("f.")
                            while pos > -1:
                                end=pos+2
                                while line[end]=='p' or line[end]=='m' or line[end]=='c':
                                    end+=1
                                off=line[pos+2:end]
                                diff=off_diff[mode][off]
                                if diff != 'c':
                                    line="%sin[i.%s%s()]%s"%(line[:pos],d,diff,line[end:])
                                else:
                                    line="%sin[i]%s"%(line[:pos],line[end:])
                                pos=line.find("f.")
                                
                            if line.find(ret) == -1:
                                print "    ",line
                            else:
                                print "    result[i]= ",line[len(ret)+line.index(ret):]
                                break
                        print """  }
  return result;
}
"""


# print at end, so functions are defined
d=dir[0]
for entry in func_db[d][:-1]:
    print "BoutReal %s(stencil &f);"%entry[0]
for d in dir:
    print "static CartesianMesh::cart_diff_lookup_table diff_lookup_%s [] = {"%d
    for entry in func_db[d]:
        print "\t{",
        t=''
        for field in entry:
            print t,field,"\t",
            t=','
        print '},'
    print "};"




for d in dir:
    warn()
    du=d.upper()
    if d != "z":
        print """const Field3D CartesianMesh::apply%sdiff(const Field3D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc ){"""%du
    else:
        print """const Field3D CartesianMesh::apply%sdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC loc ){
"""%du
    print """  cart_diff_lookup_table chosen;
  for (int i=0;;++i){
    output.write("maybe: %sp == %sp?\\n",diff_lookup_%s[i].func,func);
    if (diff_lookup_%s[i].func==func){
      chosen=diff_lookup_%s[i];
      break;
    }
    if (diff_lookup_%s[i].func==NULL){
      throw BoutException("Diff method not found!");
    }
  }
  if ( loc == CELL_%sLOW && var.getLocation() != CELL_%sLOW ){
    auto result=chosen.on(var);
    return result;
  } else if ( var.getLocation() == CELL_%sLOW  && ! (loc == CELL_CENTRE )){
    return chosen.off(var);
  } else {
    if (loc != var.getLocation() && loc != CELL_DEFAULT) throw BoutException("Didn't expect this - I wouldn't interpolate");
    auto result=chosen.norm(var);
    return result;
  }
}
"""%('%','%',d,d,d,d,du,du,du)
