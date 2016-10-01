offs={#-2:'mm',
      -1:'m',
      1:'p',
      #2:'pp'
}

#field=[2,3,'p']
field=[2,3]
dirs={}
dirs[2]=['x','y']
dirs[3]=['x','y', 'z']
dirs['p']=['x', 'z']

for f in field:
  for d in dirs[f]:
    print """
#pragma omp parallel
  for (auto i: d%s){
    d%s[i]=i.%s;
  }

"""%(f,f,d)
    for diff in offs:
        if d == 'z':
            print """
#pragma omp parallel
  for (auto i:d%s){
    auto val=d%s[i.%s%s()];
    assert(abs(val-((i.z+mesh->LocalNz*2%+d)%smesh->LocalNz))<1e-8);
  }"""%(f,f,d,offs[diff],diff,'%')
        else:
            print """
#pragma omp parallel
  for (auto i:d%s.region(RGN_NO%s)){
    auto zero=d%s[i]%+d-d%s[i.%s%s()];
    assert(abs(zero)<1e-8);
  }"""%(f,d.upper(),f,diff,f,d,offs[diff])
