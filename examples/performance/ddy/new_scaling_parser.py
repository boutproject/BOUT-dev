import csv


def read_file(filename):
    reader = csv.reader(open(filename, 'r'), delimiter='\t',
                        skipinitialspace=True)

    # Skip header
    for _, _ in zip(range(4), reader):
        continue

    from collections import OrderedDict
    case_lines = OrderedDict() #{}
    for line in reader:
        if line == []:
            break
        case_lines[line[0].rstrip('.')] = line[1]

    titles = next(reader)
    cases_weak = {col.strip(): [] for col in titles[:-1]}

    for line in reader:
        if line == []:
            break
        for title, col in zip(titles, line[:-1]):
            cases_weak[title].append(float(col))

    axis = cases_weak['Local grid']
    data = [cases_weak[x] for x in case_lines]
    labels = [case_lines[x] for x in case_lines]
    return {'axis':axis, 'data':data, 'labels':labels}

def getScan(baseName = "timing_{n}.txt", nthreads = (1,2,4,8,16,32)):
    from numpy import zeros

    dataS = []
    for n in nthreads:
        f = baseName.format(n=n)
        dataS.append(read_file(f))


    nnt = len(dataS)
    nlines = len(dataS[0]['data'])
    nval = len(dataS[0]['data'][0])
    rawDat = zeros([nnt,nval,nlines])

    for i, dat in enumerate(dataS):
        print(len(dat['data']))
        
        rawDat[i,:,:] = dat['data']

    axes = [nthreads, dataS[0]['axis']]

    return axes, rawDat, dataS[0]['labels']
