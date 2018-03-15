from __future__ import print_function

from common import *
# Define some constants, functions

import sys

# Read the stencils
stencils_raw=[]
with open("stencils_cleaned.cxx","r") as f:
    inFunc=0
    for line in f:
        if line[:5] == 'const':
            print(line)
        line=line.strip()
        # avoid to overwrite result
        line=line.replace("result","result_")
        if line == '':
            continue
        # start of function
        if line.find("{") > -1:
            if inFunc==0:
                stencils_raw.append([])
            else:
                debug("No func:",line)
            inFunc+=1
            if line.count("{")>1:
                raise RuntimeError("More than one { in \"%s\"!"%line)
        if inFunc:
            stencils_raw[-1].append(line);
        if line.find("}") > -1:
            inFunc-=1
            if line.count("}")>1:
                raise RuntimeError("More than one } in \"%s\"!"%line)

class Stencil:
    def __init__(self,body):
        self.body=body
        self.valid=True
        self.mbf="main"
        self.guards=3
        if body[0] != '':
            self.name=self.body[0].split()[1].split("(")[0]
            self.checkFlux()
            self.checkStag()
        else:
            self.name=''
            self.flux=False
            self.stag=False
        if self.name == 'DDX_CWENO3':
            debug("Skipping DDX_CWENO3 ...")
            self.valid=False
        self.setGuards()
    def setGuards(self):
        self.guards=1
        for l in self.body:
            for off in ['pp','mm','p2','m2']:
                if l.find(off)>-1:
                    self.guards=2
    def checkStag(self):
        self.stag = ( self.body[0].find("stag") > 0 )
    def checkFlux(self):
        f0=self.body[0]
        if f0.find("BoutReal") == 0 and f0.find('stencil') > 0:
            if (f0.find("BoutReal V") == -1 and f0.find("BoutReal F") == -1 ):
                # They have one stencil
                self.flux=False
            else:
                self.flux=True
        elif f0.find("&fRm") > -1:
            debug(f0,"is invalid")
            self.valid=False
        else:
            raise RuntimeError(self.body,"We did not set the flux type")


stencils=[]
for f in stencils_raw:
    sten=Stencil(f)
    if sten.valid:
        stencils.append(sten)

# Do some cleaning of the stencils bodies
for sten in stencils:
    ls=""
    for l in sten.body:
        ls+=l
    inv=ls[::-1].find("}")+1
    ls=ls[ls.find("{")+1:-inv]
    fu=ls.split(";")
    mymax=len(fu)
    i=0
    while i < mymax:
        fu[i]+=';'
        if fu[i].find("}")>-1:
            l=fu[i]
            t=l.split("}")
            for ti in t[:-1]:
                fu.insert(i,ti+"}")
                i+=1
                mymax+=1
            fu[i]=t[-1]
        i+=1
    sten.body=fu

def replace_stencil(line,sten,fname,field,mode,sten_mbf,d,update=None,z0=None):
    if update is None:
        update=sten_mbf=="main"
    pos=line.find(sten)
    part_of_offset=['p','m','c']
    for i in range(2,9):
        part_of_offset.append(str(i))
    while pos > -1:
        end=pos+len(sten)
        while line[end] in part_of_offset:
            end+=1
        off=line[pos+2:end]
        try:
            line_=line
            line=line[:pos]+get_diff(off_diff[mode][off],fname,field,d,update,z0)+line[end:]
        except:
            debug(enable=True)
            debug("Encountered error - some infos:")
            debug(line_,mode,off,sten)
            debug(off_diff)
            debug(off_diff[mode])
            raise
        pos=line.find(sten)
    return line

def parse_body(sten,field,mode,d,z0=None):
    if sten.stag == False and mode != 'norm' and sten.name:
        debug(enable=True)
        debug("Found bug - infos:")
        debug(sten.body,'name: '+sten.name,sten.stag,field,mode,d)
        assert(mode=='norm')
    body=""
    result=''
    result_=['']*2
    for line in sten.body[:-1]:
        if sten.flux:
            try:
                line=line.replace("vc","v.c");
                line=replace_stencil(line,'v.',"v_in",field,mode,sten.mbf,d,z0=z0)
                line=replace_stencil(line,'f.',"f_in",field,"norm",sten.mbf,d,z0=z0)
            except:
                debug("Some infos on error:",sten.name, sten.mbf,enable=True)
                raise
        else:
            line=replace_stencil(line,'f.',"in",field,mode,sten.mbf,d,z0=z0)
        if line.find("return") == -1:
            if sten.mbf == 'main':
                body+= "     "+line+"\n"
            else:
                toPrint=True
                for resi,res in enumerate(["result_.inner", "result_.outer"]):
                    if line.find(res) > -1:
                        tmp=line[line.index(res)+len(res):]
                        if tmp.find("=") > -1:
                            if result_[resi]!='':
                                raise RuntimeError("Did not expect another defintion of",res,
                                                   "The last one was %s = %s"%(res,result_[resi]),
                                                   "thise one is ",line)
                            result_[resi]=tmp[tmp.index("=")+1:]
                            toPrint=False
                if toPrint:
                    if line.find("=") > -1:
                        raise RuntimeError("While parsing function %s"%sten.name,"\n",
                                           sten.body,"\n",
                                           sten.mbf, "\n",
                                           "Failed to parse - unexpected line: ",line, "\n",
                                           result_, "\n",
                                           line)

        else:
            if sten.mbf == 'main':
                result=line[len("return")+line.index("return"):]+"\n"
            else:
                returned=line[len("return")+line.index("return"):]
    return [body, result, result_]

def get_for_loop_z(sten,field,stag):
    d='z'
    print('  if (Nz > 3) {')
    for d2 in dirs[field]:
        if d!=d2:
            print("  for (int %s = 0 ; %s < N%s; ++%s ){"%(d2,d2,d2,d2))
    body,result,_ = parse_body(sten,field,stag,d)
    for i in range(guards_[0]):
        print('    {')
        print('      int z=%d;'%i)
        body,result,_ = parse_body(sten,field,stag,d,z0=i+1)
        print(body)
        print("      "+get_diff('c()',"result",field,d)+"=", end=' ')
        if result:
            print(result);
        else:
            print('result_;')
        print('    }')
    print('    for (int z=%d;z<Nz-%d;++z) {'%(guards_[0],guards_[1]))
    body,result,_ = parse_body(sten,field,stag,d)
    print(body)
    print("      "+get_diff('c()',"result",field,d)+"=", end=' ')
    if result:
        print(result);
    else:
        print('result_;')
    print('    }')
    for i in range(guards_[1],0,-1):
        print('    {')
        print('      int z=Nz-%d;'%i)
        body,result,_ = parse_body(sten,field,stag,d,z0=-i)
        print(body)
        print("      "+get_diff('c()',"result",field,d)+"=", end=' ')
        if result:
            print(result);
        else:
            print('result_;')
        print('    }')
    for d2 in dirs[field]:
        if d!=d2:
            print("  }")
    print('  } else {')
    for d2 in dirs[field]:
        if d!=d2:
            print("  for (int %s = 0 ; %s < N%s; ++%s ){"%(d2,d2,d2,d2))
    print('    for (int z=0;z<Nz;++z) {')
    body,result,_ = parse_body(sten,field,stag,d,z0='secure')
    print(body)
    print("      "+get_diff('c()',"result",field,d)+"=", end=' ')
    if result:
        print(result);
    else:
        print('result_;')
    print('    }')
    for d2 in dirs[field]:
        if d!=d2:
            print("  }")
    print('  }')

def get_for_loop(d,mode,field,guards,sten_name ):
    if sten_name == "main":
        print('#if CHECK > 0')
        print('  if (msh->%sstart < %d){'%(d,max(guards)))
        print('    throw BoutException("Cannot compute derivative - need at least %d guard cells in %s direction!");'%(max(guards),d.upper()))
        print('  }')
        print('#endif')
        dp=guards[0]
        dm=-guards[1]
        for d2 in dirs[field]:
            if d==d2:
                print("  for (int %s = %d; %s < N%s%+d; ++%s ){"%(d,dp,d,d,dm,d))
            else:
                print("  for (int %s = 0 ; %s < N%s; ++%s ){"%(d2,d2,d2,d2))
    else:
        if sten_name == 'forward':
            print("    int "+d,"=%d ;"%(guards_[0]-1))
        else:
            print("    "+d,"=N%s"%d,"-%d ;"%(guards_[1]))
            if guards_[1]>2:
                raise RuntimeError(guards_,"To many guards")
        for d2 in perp_dir[field][d]:
            print("    for (int "+d2,"=0; "+d2,"< N"+d2,";++"+d2,") {")


def get_for_end(d,field, sten_name):
    if sten_name == 'main':
        for d2 in dirs[field]:
            print("  }")
    else:
        for d2 in perp_dir[field][d]:
            print("  }")

def get_diff(diff,fname,field,d,update=False,z0=None):
    global use_field_operator
    if use_field_operator:
        ret=fname+"("
    else:
        ret=fname+"_ptr["
    for d2 in dirs[field][1:]:
        if not use_field_operator:
            ret+="("
    first=True
    for d2 in dirs[field]:
        if not first:
            if use_field_operator:
                ret+=','
            else:
                ret+=')*N%s + '%d2
        first=False
        if (d2 == d):
            #do ugly stuff
            global guards_
            diffd=diff2[diff]
            if update:
                if abs(diffd) > 2:
                    raise RuntimeError("We do not expect to have more then 2 guard cells")
                if diffd > guards_[1]:
                    guards_[1]=diffd
                elif -diffd > guards_[0]:
                    guards_[0]=-diffd
                #else:
                    #print diffd,guards_;
            if d=='z' and z0 is not None:
                # We want to generate code for the interp_to case with wrapping
                if z0 == 'secure':
                    ret+="+(("+d+"%+d"%(diffd)
                    if diffd < 0:
                        ret+='+%d*Nz'%(-diffd)
                    ret+=')%Nz)'
                elif z0 > 0 and -diffd >= z0:
                    ret+=d2+"%+d+Nz"%(diffd)
                elif z0 < 0 and -diffd <= z0:
                    ret+=d2+"%+d-Nz"%(diffd)
                else:
                    ret+=d2+"%+d"%(diffd)
            else:
                ret+=d2+"%+d"%(diffd)
        else:
            ret+=d2
    if use_field_operator:
        ret+=")\n    \t\t\t\t"
    else:
        ret+="]\n    \t\t\t\t"
    return ret

def get_pointer(field, field_type,const):
    if const:
        print('  checkData(%s);'%field)
        print("  const BoutReal * __restrict__", end=' ')
    else:
        print("  BoutReal * __restrict__", end=' ')
    print("%s_ptr = &"%field,field,"(", end=' ')
    first=True
    for d in dirs[field_type]:
        if not first:
            print(',', end=' ')
        first=False
        print("0", end=' ')
    print(");")


def gen_functions_normal(to_gen):
    for ftg in to_gen:
        mode=ftg.stag_mode
        d=ftg.d
        field=ftg.field
        name=ftg.name
        flux=ftg.flux
        warn()
        print("static void "+ftg.name+"_"+ftg.field.lower()+"(BoutReal * __restrict__ result_ptr,", end=' ')
        if flux:
            print("const BoutReal * __restrict__ v_in_ptr,", end=' ')
            print("const BoutReal * __restrict__ f_in_ptr,", end=' ')
        else:
            print("const BoutReal * __restrict__ in_ptr,", end=' ')
        print("Mesh * msh) {")
        for d2 in dirs[field]:
            print("  const int N%s = msh->LocalN%s;"%(d2,d2))
        for func in stencils:
            if func.name==ftg[4]:
                ftg.setSten(func)
                assert(ftg.sten is not None)
        # only reset before main
        global guards_
        guards_=[0,0]
        if ftg.sten is None:
            debug("infos on error:",enable=True)
            debug(ftg,ftg[:])
            for func in stencils:
                debug(func.name)
            debug(ftg.sten)
            debug(ftg)
            raise RuntimeError("error unexpected: sten is None for sten_name %s !"%sten_name)
        guards=ftg.sten.guards
        if d=='z':
            get_for_loop_z(ftg.sten,field,mode)
        else:
            body, result, result_ = parse_body(ftg.sten,field,mode,d)
            get_for_loop(d,mode,field,guards_,"main")
            if body+result == '':
                raise RuntimeError(ftg.sten.body)
            print(body)
            print("      "+get_diff('c()',"result",field,d)+"=", end=' ')
            if result:
                print(result);
            else:
                print('result_;')
            get_for_end(d,field,"main")
            
        print("}")
        warn()
        print("static",field,name,"(const",field, end=' ')
        if flux:
            print("&v_in, const",field,"&f_in){")
            print("  Mesh * msh = v_in.getMesh();")
            print("  ASSERT1(msh == f_in.getMesh());")
        else:
            print("&in){")
            print("  Mesh * msh = in.getMesh();")
        print('  output_debug.write("Using method %s!\\n");'%name)
        if d=='z':
            print('  if (msh->LocalN%s == 1) {'%(d))
            print('    return %s(0.,msh);'%field)
            print('  }')
        print('#if CHECK > 0')
        print('  if (msh->LocalN%s < %d) {'%(d,sum(guards_)+1))
        print('    if (msh->xstart == 0) {')
        print('      //return %s(0.,msh);'%field)
        print('      %s result{msh};'%field)
        print('      result=0;')
        print('      return result;')
        print('    } else {')
        print('      throw BoutException("AiolosMesh::%s - Not enough guards cells to take derivative!");'%(name))
        print('    }')
        print('  }')
        print('#endif')
        print(" ",field,"result(msh);")
        print("  result.allocate();")
        get_pointer("result",field,False)
        if flux:
            get_pointer("v_in",field,True)
            get_pointer("f_in",field,True)
        else:
            get_pointer("in",field,True)
        if flux:
            print("  "+name+"_"+field.lower()+"(result_ptr,v_in_ptr,f_in_ptr, msh);")
        else:
            print("  "+name+"_"+field.lower()+"(result_ptr,in_ptr, msh);")
        if mode == "on":
            print("  result.setLocation(CELL_%sLOW);"%d.upper())
        elif mode == "off":
            print("  result.setLocation(CELL_CENTRE);")
        else:
            if flux:
                print("  result.setLocation(f_in.getLocation());")
            else:
                print("  result.setLocation(in.getLocation());")
        print("""  checkData(result);
  return result;
}
""")


field=fields[0]#3d
use_field_operator=True
def get_interp_vals(order,pos):
    import numpy as np
    rhs=np.zeros(order)
    rhs[0]=1
    mat=np.zeros((order,order))
    x0=-pos+(order-1.)/2.
    for i in range(order):
        x=x0-i;
        for j in range(order):
            mat[j,i]=x**j*np.math.factorial(j)
    #debug(mat)
    facs=np.dot(np.linalg.inv(mat),rhs)
    return facs
def get_interp_sten(order,pos):
    if order%2:
        print("Cannot handle uneven order!", file=sys.stderr)
        exit(4)
    oh=order/2
    if pos:
        oh-=pos/abs(pos)
    vals=get_interp_vals(order,pos)
    ret=""
    first=True
    for i in range(order):
        ret+="%+.5e*"%vals[i]
        #debug(ret)
        if i < oh:
            ret+='f.m'+(str(int(oh-i)) if oh-i > 1 else "")
        else:
            ret+='f.p'+(str(int(i-oh+1)) if i-oh+1 > 1 else "" )
    #debug(vals)
    #debug(ret)
    return ret+" ;"


interp=['',"return "+get_interp_sten(4,0),'']
# Hard coded version:
#["return ( 9.*(s.m + s.p) - s.mm - s.pp ) / 16.;"]
for mode in ['on','off']:
  for order in [4]:
    for d in dirs[field]:
        global guards_
        guards_=[0,0]
        line=interp[1]
        sten_name="main"
        line=replace_stencil(line,'f.',"in",field,mode,sten_name,d)
        print("static void interp_to_%s_%s_%s("%(mode,field,d), end=' ')
        if use_field_operator:
            print(field+"& result, const "+field+" & in,", end=' ')
        else:
            print("BoutReal * __restrict__ result_ptr, const BoutReal * __restrict__ in_ptr,", end=' ')
        print(" Mesh * msh ){")
        for d2 in dirs[field]:
            print("  const int N%s = msh->LocalN%s;"%(d2,d2))
        if d == 'z':
            sten=Stencil(['',interp[1],''])
            get_for_loop_z(sten,field,mode)
        else:
            body= "    "+get_diff('c()',"result",field,d,update=True)+"= "+line[len("return")+line.index("return"):]+"\n"
            get_for_loop(d,mode,field,guards_,sten_name)
            print(body)
            guards__=guards_
            get_for_end(d,field,sten_name)
            sten_names=["forward","backward"]
            for sten_name in sten_names:
                sten_name_index=sten_names.index(sten_name)
                sign=-1
                if sten_name == "forward":
                    sign=1
                _sign=sign
                if d=='z':
                    sign=0
                get_for_loop(d,mode,field,guards_,sten_name)
                print("      "+get_diff('c()',"result",field,d)+"=", end=' ')
                print(replace_stencil(get_interp_sten(4,sign),'f.',"in",field,mode,sten_name,d,False,z0=((guards_[sten_name_index])*_sign)))
                guards_=guards__
                if order/2 > 1:
                    if (sten_name == 'backward' and mode == 'on') or \
                       (sten_name == 'forward' and mode == 'off'):
                        pass# dont do anything ...
                    else:
                        if sten_name == 'forward':
                            print("        "+get_diff('m()',"result",field,d)+"=", end=' ')
                        else:
                            print("        "+get_diff('p()',"result",field,d)+"=", end=' ')
                        print(replace_stencil(get_interp_sten(4,sign*2),'f.',"in",field,mode,sten_name,d,False,z0=((guards_[sten_name_index])*_sign)))
                        guards_=guards__
                get_for_end(d,field,sten_name)
        print("}")

print("const Field3D AiolosMesh::interp_to_do(const Field3D &f, CELL_LOC loc) const {")
print("  Mesh * msh = f.getMesh();")
print("  Field3D result(msh);")
print("  result.allocate();")
print("  if (f.getLocation() != CELL_CENTRE){")
print("    // we first need to go back to centre before we can go anywhere else")
print("    switch (f.getLocation()){")
for d in dirs[field]:
    print("    case CELL_%sLOW:"%d.upper())
    if use_field_operator:
        print("      interp_to_off_%s_%s(result,f,msh);"%(field,d))
    else:
        print("      interp_to_off_%s_%s(&result(0,0,0),&f(0,0,0),msh);"%(field,d))
    print("      result.setLocation(CELL_CENTRE);")
    print("      // return or interpolate again")
    print("      return interp_to(result,loc);")
    print("      break;")
print("    default:")
print('      throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",strLocation(loc));')
print("    }")
print("  }")
print("  // we are in CELL_CENTRE and need to go somewhere else ...")
print("  switch (loc){")
for d in dirs[field]:
    print("    case CELL_%sLOW:"%d.upper())
    if use_field_operator:
        print("      interp_to_on_%s_%s(result,f,msh);"%(field,d))
    else:
        print("      interp_to_on_%s_%s(&result(0,0,0),&f(0,0,0),msh);"%(field,d))
    print("      result.setLocation(CELL_%sLOW);"%d.upper())
    print("      // return or interpolate again")
    print("      return interp_to(result,loc);")
    print("      break;")
print("    default:")
print('      throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",strLocation(loc));')
print("    }")
print("}")

use_field_operator=False
