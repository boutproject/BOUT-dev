import csv


def read_file(filename):
    reader = csv.reader(open(filename, 'r'), delimiter='\t',
                        skipinitialspace=True)

    # Skip header
    for _, _ in zip(range(4), reader):
        continue

    case_lines = {}
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

    return case_lines, cases_weak, titles
