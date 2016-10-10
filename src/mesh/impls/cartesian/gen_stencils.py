from common import *
# Define some constants, functions

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
                
# for f in functions:
#     if f.flux == False:
#             if f.stag:
#                 modes=['on','off']
#             else:
#                 modes=['norm']
#             for d in dir:
#                 tmp="cart_diff_%s_%s_%%s"%(d,f.name)
#                 for mode in modes:
#                         warn()
#                         print """Field3D %s(const Field3D in){
#   Field3D result;
#   //result.allocate();
#   result=0;"""%(tmp%mode)
#                         if d != 'z':
#                             print """  for (auto i: result.region(RGN_NO%s)){"""%d.upper()
#                         else:
#                             print """  for (auto i: result){"""
#                         for line in f.body[1:]:
#                             ret="return"
#                             pos=line.find("f.")
#                             while pos > -1:
#                                 end=pos+2
#                                 while line[end]=='p' or line[end]=='m' or line[end]=='c':
#                                     end+=1
#                                 off=line[pos+2:end]
#                                 diff=off_diff[mode][off]
#                                 if diff != 'c':
#                                     line="%sin[i.%s%s()]%s"%(line[:pos],d,diff,line[end:])
#                                 else:
#                                     line="%sin[i]%s"%(line[:pos],line[end:])
#                                 pos=line.find("f.")
                                
#                             if line.find(ret) == -1:
#                                 print "    ",line
#                             else:
#                                 print "    result[i]= ",line[len(ret)+line.index(ret):]
#                                 break
#                         print """  }
#   return result;
# }
# """

                
def gen_functions_normal(to_gen):
    for f_ar in to_gen:
        mode=f_ar[3]
        #for d in dir:
        d=f_ar[2]
        field=f_ar[1]
        #tmp="cart_diff_%s_%s_%%s"%(d,f.name)
        #for mode in modes:
        name=f_ar[0]
        warn()
        print "const",field,name,"(const",field,"&in){"
        print " ",field,"result;"
        print "  result.allocate();"
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
                if d != 'z':
                    print """  for (auto i: result.region(RGN_NO%s)){"""%d.upper()
                else:
                    print """  for (auto i: result){"""
            elif sten_name!='main':
                if sten_name == 'forward':
                    print "  if (mesh->%sstart > 0){"%d
                    print "    DataIterator i(0,mesh->LocalNx,0,mesh->LocalNy,0,mesh->LocalNz);"
                    print "    i."+d,"=mesh->%sstart"%d,"-1 ;"
                if sten_name == 'backward':
                    print "    i."+d,"=mesh->%send"%d,"+1 ;"
                for d2 in perp_dir[field][d]:
                    print "    for (i."+d2,"=0; i."+d2,"< mesh->LocalN"+d2,";++i."+d2,") {"
            sten=stencils[sten_name]
            if sten is None:
                print f_ar
                for func in functions:
                    print func.name
                print "#error unexpected: sten is None!"
                exit(1)
            for line in sten.body[1:]:
                ret="return"
                pos=line.find("f.")
                while pos > -1:
                    end=pos+2
                    part_of_offset=['p','m','c']
                    for i in range(2,9):
                        part_of_offset.append(str(i))
                    while line[end] in part_of_offset:
                        end+=1
                    off=line[pos+2:end]
                    diff=off_diff[mode][sten_name][off]
                    if diff != 'c':
                        line=line[:pos]+"in[i.%s%s()]"%(d,diff)+line[end:]
                    else:
                        line="%sin[i]%s"%(line[:pos],line[end:])
                    pos=line.find("f.")

                if line.find(ret) == -1:
                    print "     ",line
                else:
                    if sten_name == 'main':
                        print "    result[i]= ",line[len(ret)+line.index(ret):]
                    else:
                        returned=line[len(ret)+line.index(ret):]
                    break
            if sten_name != 'main':
                print "      result[i]=result_.inner;"
                print "      if (mesh->%sstart >1 ){"%d
                if sten_name == 'forward':
                    print "        result[i.%sm()]=result_.outer;"%d
                else:
                    print "        result[i.%sp()]=result_.outer;"%d
                print "      }"
                for d2 in perp_dir[field][d]:
                    print "    }"
            else:
                print "  }"
            if sten_name == 'backward':
                print "  }"
            
        print """  return result;
}
"""


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
