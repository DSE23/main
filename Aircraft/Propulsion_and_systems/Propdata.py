"""
Name: Propeller Data Reader
Department: Propulsion and Aircraft Systems
Last updated: 19/06/2018 10:54 by Ties
"""

"""
This file provides an easy way to extract data from the DataReal.txt file
"""

import csv

with open('DataReal.txt', 'r') as Data:
    reader = csv.DictReader(Data, delimiter=" ")
    data = []

    for row in reader:
        if row != 0:
            data.append(row)

    print(data[:][3])