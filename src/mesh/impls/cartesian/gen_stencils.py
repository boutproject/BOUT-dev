from common import *
# Define some constants, functions

import sys

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
    guards=3

functions_=functions
functions=[]
for f in functions_:
    func=function()
    func.body=f
    func.name=f[0].split()[1].split("(")[0]
    if f[0].find("stag") > 0:
        # They are staggered
        func.stag=True
    else:
        func.stag=False
    if f[0].find("BoutReal") == 0 and f[0].find('stencil') > 0:
        #for i in f:
        #    print i
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
        functions.append(func)
    if f[0].find("Mesh::boundary_derivs_pair") == 0 and f[0].find('stencil') > 0:
        if (f[0].find("Mesh::boundary_derivs_pair V") == -1 and f[0].find("Mesh::boundary_derivs_pair F") == -1 ):
            # not a flux functions - has only on stencil as argument
            func.flux=False
        else:
            func.flux=True
        functions.append(func)
null_func=function()
null_func.name='NULL'
null_func.body=["NULL","result_.inner = 0;","result_.outer = 0;","}"]
functions.append(null_func)
for f in functions:
    f.guards=1
    for l in f.body:
        for off in ['pp','mm','p2','m2']:
            if l.find(off)>-1:
                f.guards=2
                #print >> sys.stderr , f.name ,"has 2 guard cells"
    #if f.guards == 1:
        #print >> sys.stderr , f.name ,"has 1 guard cells"
        
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


def replace_stencil(line,sten,fname,field,mode,sten_name,off,d):
    pos=line.find(sten)
    while pos > -1:
        end=pos+len(sten)
        part_of_offset=['p','m','c']
        for i in range(2,9):
            part_of_offset.append(str(i))
        while line[end] in part_of_offset:
            end+=1
        off=line[pos+2:end]
        try:
            line_=line
            line=line[:pos]+get_diff(off_diff[mode][sten_name][off],fname,field,d)+line[end:]
        except:
            print >>sys.stderr,line_
            raise
        pos=line.find(sten)
    return line


def get_for_loop(d,mode,field,guards,sten_name ):
    if sten_name == 'main':
        if d != 'z':
            dp=guards
            dm=-guards
            if mode == 'on':
                dm=-(guards-1)
            elif mode == 'off':
                dp=(guards-1)
            # if 'x' in dirs[field]:
            #     xmax='mesh->LocalNx%+d'%(-1+dxm)
            # else:
            #     xmax="0"
            # if 'y' in dirs[field]:
            #     ymax='mesh->LocalNy%+d'%(-1+dym)
            # else:
            #     ymax="0"
            # if 'z' in dirs[field]:
            #     zmax='mesh->LocalNz-1'
            # else:
            #     zmax="0"
            #print "  for (DataIterator i(%s, %s,"%(dxp,xmax),
            #print "%s, %s, 0, %s)"%(dyp,ymax,zmax),
            #print "; !i.done() ; ++i) {"
        if d=='z':
            d='ignore'
        for d2 in dirs[field]:
            if d==d2:
                print "  for (int %s = %d; %s < N%s%+d; ++%s ){"%(d,dp,d,d,dm,d)
            else:
                print "  for (int %s = 0 ; %s < N%s; ++%s ){"%(d2,d2,d2,d2)
        #else:
        #    print """  for (auto i: result){"""
    else:
        if sten_name == 'forward':
            print "  if (mesh->%sstart > 0){"%d
            #print "    DataIterator i(0,mesh->LocalNx,0,mesh->LocalNy,0,mesh->LocalNz);"
        bp=guards
        if ( sten_name == 'backward' and mode == 'on' ) or \
           ( sten_name == 'forward' and mode == 'off' ):
            if guards > 1:
                bp=guards-1;
        if sten_name == 'forward':
            print "    int "+d,"=%d ;"%(bp-1)
        else:
            print "    "+d,"=mesh->LocalN%s"%d,"-%d ;"%bp
        for d2 in perp_dir[field][d]:
            print "    for (int "+d2,"=0; "+d2,"< mesh->LocalN"+d2,";++"+d2,") {"
            

def get_for_end(d,field, sten_name):
    if sten_name != 'main':
        for d2 in perp_dir[field][d]:
            print "    }"
    else:
        for d2 in dirs[field]:
            print "  }"
    if sten_name == 'backward':
        print "  }"

def get_diff(diff,fname,field,d):
    #diff=
    #print >> sys.stderr, diff
    ret=fname+"_ptr["
    for d2 in dirs[field][1:]:
        ret+="("
    first=True
    for d2 in dirs[field]:
        if not first:
            ret+=')*N%s + '%d2
        first=False
        if (d2 == d):
            #do smart stuff
            ret+=d2+"%+d"%(diff2[diff])
        else:
            ret+=d2
    ret+="]"
    return ret
    # if diff != 'c()':
    #     return "%s[i.%s%s]"%(fname,d,diff)
    # else:
    #     return "%s[i]"%(fname)

def get_pointer(field, field_type,const):
    if const:
        print "  const BoutReal * __restrict__",
    else:
        print "  BoutReal * __restrict__",
    print "%s_ptr = &"%field,field,"(",
    first=True
    for d in dirs[field_type]:
        if not first:
            print ',',
        first=False
        print "0",
    print ");";
        

def gen_functions_normal(to_gen):
    import sys
    #for f_ar in to_gen:
    #    print >>sys.stderr,f_ar[0],f_ar[1]
    #exit(1)
    for f_ar in to_gen:
        mode=f_ar[3]
        #for d in dir:
        d=f_ar[2]
        field=f_ar[1]
        #tmp="cart_diff_%s_%s_%%s"%(d,f.name)
        #for mode in modes:
        name=f_ar[0]
        flux=f_ar[-1]
        warn()
        print "static void "+name+"_"+field.lower()+"(BoutReal * __restrict__ result_ptr,",
        if flux:
            print "const BoutReal * __restrict__ v_in_ptr,",
            print "const BoutReal * __restrict__ f_in_ptr) {"
        else:
            print "const BoutReal * __restrict__ in_ptr) {"
        for d2 in dirs[field]:
            print "  const int N%s = mesh->LocalN%s;"%(d2,d2)
        stencils={'main':None,
                  'forward':None,
                  'backward':None}
        for func in functions:
            if func.name==f_ar[4]:
                stencils['main']=func
            if func.name==f_ar[5]:
                stencils['forward']=func
            if func.name==f_ar[6]:
                stencils['backward']=func
        #start with main file:
        todo=['main','forward','backward']
        if d=='z':
            todo=[todo[0]]
        for sten_name in todo:
            if sten_name=='main':
                guards=stencils[sten_name].guards #numGuards[f_ar[4]]
            if sten_name=='backward' and mode=='on' and guards ==1:
                print "  }"
                continue;
            if sten_name=='forward' and mode=='off' and guards ==1:
                print "  if (mesh->%sstart > 0){"%d
                #print "    DataIterator i(0,mesh->LocalNx,0,mesh->LocalNy,0,mesh->LocalNz);"
                print "    int",d,";"
                continue;
            get_for_loop(d,mode,field,guards,sten_name)
            sten=stencils[sten_name]
            if sten is None:
                import sys
                print >>sys.stderr,f_ar
                for func in functions:
                    print >>sys.stderr,func.name
                print >>sys.stderr,stencils
                print >>sys.stderr,"#error unexpected: sten is None for sten_name %s !"%sten_name
                exit(1)
            result_=['',''] # for foward/backward
            for line in sten.body[1:-1]:
                if flux:
                    line=replace_stencil(line,'v.',"v_in",field,mode,sten_name,off,d)
                    line=replace_stencil(line,'f.',"f_in",field,"norm",sten_name,off,d)
                else:
                    line=replace_stencil(line,'f.',"in",field,mode,sten_name,off,d)

                if line.find("return") == -1:
                    if sten_name == 'main':
                        print "     ",line
                    else:
                        toPrint=True
                        resl=["result_.inner", "result_.outer"]
                        for resi in range(2):
                            res=resl[resi]
                            if line.find(res) > -1:
                                tmp=line[line.index(res)+len(res):]
                                if tmp.find("=") > -1:
                                    if result_[resi]!='':
                                        import sys
                                        print >> sys.stderr ,"Did not expect another defintion of",res
                                        print >> sys.stderr ,"The last one was %s = %s"%(res,result_[resi])
                                        print >> sys.stderr ,"thise one is ",line
                                        exit(1)
                                    result_[resi]=tmp[tmp.index("=")+1:]
                                    toPrint=False
                        if toPrint:
                            if line.find("=") > -1:
                                import sys
                                print >> sys.stderr ,"Failed to parse - unexpected line: ",line
                                print >> sys.stderr ,result_
                                exit(1)
                                    
                else:
                    if sten_name == 'main':
                        print "    "+get_diff('c()',"result",field,d)+"= ",line[len("return")+line.index("return"):]
                    else:
                        returned=line[len("return")+line.index("return"):]
                    #break
            if sten_name != 'main':
                if result_[0] != '':
                    print "      "+get_diff('c()',"result",field,d)+"="+result_[0]
                else:
                    import sys
                    print "      "+get_diff('c()',"result",field,d)+"=result_.inner;"
                #print "      if (mesh->%sstart >1 ){"%d
                
                if guards > 1:
                    #print >> sys.stderr, stencils[todo[0]].body[0]
                    #print >> sys.stderr, stencils[todo[0]].guards
                    if ( sten_name == 'backward' and mode == 'on' ) or \
                       ( sten_name == 'forward' and mode == 'off' ):
                        pass# dont do anything ...
                    else:
                        if sten_name == 'forward':
                            print "        "+get_diff('m()',"result",field,d)+"=" ,
                        else:
                            print "        "+get_diff('p()',"result",field,d)+"=" ,
                        if result_[1] != '':
                            print result_[1]
                        else:
                            print "result_.outer;"
                #print "      }"
            get_for_end(d,field,sten_name)
            
        print "}"
        warn()
        print "static",field,name,"(const",field,
        if flux:
            print "&v_in, const",field,"&f_in){"
        else:
            print "&in){"
        print ' output.write("Using method %s!\\n");'%name
        print " ",field,"result;"
        print "  result.allocate();"
        get_pointer("result",field,False)
        if flux:
            get_pointer("v_in",field,True)
            get_pointer("f_in",field,True)
        else:
            get_pointer("in",field,True)
        if flux:
            print "  "+name+"_"+field.lower()+"(result_ptr,v_in_ptr,f_in_ptr);"
        else:
            print "  "+name+"_"+field.lower()+"(result_ptr,in_ptr);"
        print """  return result;
}
"""
    #exit(1)

# # print at end, so functions are defined
# d=dir[0]
# for entry in func_db[d][:-1]:
#     print "BoutReal %s(stencil &f);"%entry[0]
# for d in dir:
#     print "static CartesianMesh::cart_diff_lookup_table diff_lookup_%s [] = {"%d
#     for entry in func_db[d]:
#         print "\t{",
#         t=''
#         for field in entry:
#             print t,field,"\t",
#             t=','
#         print '},'
#     print "};"




# for d in dir:
#     warn()
#     du=d.upper()
#     if d != "z":
#         print "const Field3D CartesianMesh::apply%sdiff(const Field3D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc ){"%du
#     else:
#         print "const Field3D CartesianMesh::apply%sdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC loc ){"%du
#     print """  cart_diff_lookup_table chosen;
#   for (int i=0;;++i){
#     output.write("maybe: %sp == %sp?\\n",diff_lookup_%s[i].func,func);
#     if (diff_lookup_%s[i].func==func){
#       chosen=diff_lookup_%s[i];
#       break;
#     }
#     if (diff_lookup_%s[i].func==NULL){
#       throw BoutException("Diff method not found!");
#     }
#   }
#   if ( loc == CELL_%sLOW && var.getLocation() != CELL_%sLOW ){
#     auto result=chosen.on(var);
#     return result;
#   } else if ( var.getLocation() == CELL_%sLOW  && ! (loc == CELL_CENTRE )){
#     return chosen.off(var);
#   } else {
#     if (loc != var.getLocation() && loc != CELL_DEFAULT) throw BoutException("Didn't expect this - I wouldn't interpolate");
#     auto result=chosen.norm(var);
#     return result;
#   }
# }
# """%('%','%',d,d,d,d,du,du,du)
