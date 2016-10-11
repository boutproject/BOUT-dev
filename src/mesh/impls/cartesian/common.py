# Define some constants, functions
def warn():
    print "// This file is auto-generated - do not edit!"

fields=['Field3D','Field2D']
dirs=dict()
dirs['Field3D']=['x', 'y', 'z']
dirs['Field2D']=['x', 'y']
perp_dir=dict()
for field in fields:
    perp_dir[field]=dict()
    for d in dirs[field]:
        perp_dir[field][d]=[]
        for d2 in dirs[field]:
            if d != d2:
                perp_dir[field][d].append(d2)

off_diff=dict()
for stag in ['on','off','norm']:
    off_diff[stag]=dict()
    for mfb in ['main','forward','backward']:
        curr=dict()
        if mfb=='main':
            if stag=='on':
                curr={'ppp':'ppp',
                      'pp' :'pp',
                      'p'  :'p',
                      'm'  :'c',
                      'mm' :'m',
                      'mmm':'mm'}
            if stag=='off':
                curr={'ppp':'pp',
                      'pp' :'p',
                      'p'  :'c',
                      'm'  :'m',
                      'mm' :'mm',
                      'mmm':'mmm'}
            if stag=='norm':
                curr={'ppp':'ppp',
                      'pp' :'pp',
                      'p'  :'p',
                      'c'  :'c',
                      'm'  :'m',
                      'mm' :'mm',
                      'mmm':'mmm'}
            for d in curr:
                curr[d]=curr[d]+'()'
                if curr[d]=='mmm()':
                    curr[d]='m(3)'
                if curr[d]=='ppp()':
                    curr[d]='p(3)'
        else:
            mylist=[]
            mylist2=[]
            for i in range(7,1,-1):
                mylist.append('m%d'%i)
                mylist2.append('m(%d)'%i)
            mylist.append('m')
            mylist.append('c')
            mylist.append('p')
            mylist2.append('m()')
            mylist2.append('c()')
            mylist2.append('p()')
            for i in range(2,8):
                mylist.append('p%d'%i)
                mylist2.append('p(%d)'%i)
            if stag == 'norm' or \
               (mfb=='backward' and stag=='off' ) or \
               (mfb=='forward' and stag=='on'   ):
                for i in range(len(mylist)):
                    curr[mylist[i]]=mylist2[i];
            elif mfb == 'forward' and stag=='off':
                for i in range(len(mylist)-1):
                    curr[mylist[i]]=mylist2[i+1]
            elif mfb == 'backward' and stag=='on':
                for i in range(len(mylist)-1):
                    curr[mylist[i+1]]=mylist2[i]
        off_diff[stag][mfb]=curr

