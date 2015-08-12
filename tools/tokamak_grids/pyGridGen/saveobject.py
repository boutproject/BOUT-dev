from future import standard_library
standard_library.install_aliases()
import pickle as pickle


def saveobject(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


