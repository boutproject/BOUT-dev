from random import randrange as rand
class Field:
    n=''
    d=[]
    i=''
    def __init__(self,name,dirs,idn):
        self.n=name
        self.d=dirs
        self.i=idn

    def get(self,name):
        if self.i=='real':
            return name
        elif self.i=='f2d':
            return '%s[y+x*ny]'%name
        else:
            return '%s[z+nz*(y+ny*x)]'%(name)

    def posibles(self):
        ret=["%g"%(rand(1)+1)]
        for d1 in self.d:
            ret.append("%s+%g"%(d1,1+rand(1)))
            for d2 in self.d:
                if d1 < d2:
                    ret.append("%s*%s+%g"%(d1,d2,rand(1)+1))
                    for d3 in self.d:
                        if d2 < d3:
                            ret.append("%s*%s*%s+%g"%(d1,d2,d3,rand(1)+1))
        return ret
            



f3d =Field('Field3D' ,['x','y','z'],'f3d')
f2d =Field('Field2D' ,['x','y'    ],'f2d')
real=Field('BoutReal',[           ],'real')
fields=[f3d,f2d,real]

def mymax(f1,f2):
    if f1.i==f2.i:
        return f1
    elif f1.i == 'real':
        return f2
    elif f2.i == 'real':
        return f1
    else:
        return f3d
def mymin(f1,f2):
    if (len(f1.d) < len(f2.d)):
        return f1
    else:
        return f2

ops='*/+-'
op_names={
    '*' : "Multiplication",
    '/' : "Division",
    '+' : "Addition",
    '-' : "Substraction"
}

for field in fields:
    print "%s %s(const string &value){"%(field.n,field.i)
    if field.i == 'f3d':
        print "  return FieldFactory::get()->create3D(value,NULL,mesh,CELL_CENTRE,0);"
    elif field.i == 'f2d':
        print "  return FieldFactory::get()->create2D(value,NULL,mesh,CELL_CENTRE,0);"
    elif field.i == 'real':
        print "  return atof(value.c_str());"
    print "}"
    print

for f1 in fields:
    for f2 in fields:
        if f1.i == f2.i == 'real':
            continue
        out=mymax(f1,f2)
        for op in ops:
            print "void test_%s_%s_%s(){"%(f1.i,f2.i,op_names[op])
            for i1 in f1.posibles():
                for i2 in f2.posibles():
                    print "  compare(%s(\"(%s)%s(%s)\")"%(out.i,i1,op,i2)
                    print '  \t,%s("%s")%s%s("%s")'%(f1.i,i1,op,f2.i,i2)
                    print '  \t,"%s of %s(%s) and %s(%s)");'%(op_names[op],f1.n,i1,f2.n,i2)
                    #print '}'
                    print
            print "}"
            print

with open("test_fieldops.cxx","w") as out:
    for f1 in fields:
        for f2 in fields:
            if f1.i == f2.i == 'real':
                continue
            for op in ops:
                print >>out ,"  test_%s_%s_%s();"%(f1.i,f2.i,op_names[op])
