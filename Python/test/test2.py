import sys

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path

from pickler import storeData, loadData  # Imports pickler functions

data = {"list" : [1, 2, 3]}

storeData(data, "db")
loadData("db", showKeys=True)
