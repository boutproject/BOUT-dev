# Define some constants, functions
def warn():
    print "\n// This file is auto-generated - do not edit!"

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
            if stag=='off':
                curr={'ppp':'ppp',
                      'pp' :'pp',
                      'p'  :'p',
                      'm'  :'c',
                      'mm' :'m',
                      'mmm':'mm'}
            if stag=='on':
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

diff2={'mm()': -2,
       'm()' : -1,
       'c()' :  0,
       'p()' :  1,
       'pp()':  2}
for i in range(10):
    diff2["m(%d)"%i]=-i
    diff2["p(%d)"%i]=i
numGuards={ 'DDX_C2'       : 1,
            'DDX_C4'       : 2,
            'DDX_U1'       : 1,
            'DDX_U2'       : 2,
            'DDX_U4'       : 2,
            'DDX_CWENO2'       : 1,
            'DDX_WENO3'       : 2,
            'DDX_S2'       : 2,
            'DDX_NND'      : 2,
            'DDX_C2_stag'  : 1,
            'DDX_C4_stag'  : 2,
            'DDX_U1_stag'  : 1,
            'DDX_U2_stag'  : 2,
            'DDX_U4_stag'  : 2,
            'DDX_WENO2_stag'  : 1,
            'DDX_WENO3_stag'  : 2,
            'DDX_S2_stag'  : 2,
            'DDX_NND_stag' : 2,
            'D2DX2_C2'       : 1,
            'D2DX2_C4'       : 2,
            'D2DX2_U1'       : 1,
            'D2DX2_U2'       : 2,
            'D2DX2_U4'       : 2,
            'D2DX2_WENO2'       : 1,
            'D2DX2_WENO3'       : 2,
            'D2DX2_S2'       : 2,
            'D2DX2_NND'      : 2,
            'D2DX2_C2_stag'  : 1,
            'D2DX2_C4_stag'  : 2,
            'D2DX2_U1_stag'  : 1,
            'D2DX2_U2_stag'  : 2,
            'D2DX2_U4_stag'  : 2,
            'D2DX2_WENO2_stag'  : 1,
            'D2DX2_WENO3_stag'  : 2,
            'D2DX2_S2_stag'  : 2,
            'D2DX2_NND_stag' : 2,
            }
