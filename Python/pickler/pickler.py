import pickle
import os
from datetime import datetime
import shutil

#################################
# --- Data storage function --- #
#################################

def storeData(data : dict, fileName : str, cwd=os.getcwd()):
    dest = cwd + "/" + fileName  # Destination of data file
    print("Storing file at " + dest)

    # Backups old db if such exists and then overwrites the old one after if-clause
    if os.path.isfile(dest):
        src = dest  # For clarity renames 'dest', the path of already existing file, to 'src'
        
        backupName = fileName + datetime.today().strftime('_%Y-%m-%d_%H:%M:%S')  # New name of file
        backupDir = cwd + "/backups/"  # Directory of backups
        os.makedirs(os.path.dirname(backupDir), exist_ok=True)  # Checks whether said directory exists. Otherwise it is created.
        newDest = backupDir + backupName  # Destination of backup file
        shutil.copy(src, newDest)  # Creates backup of file

        print(f"File already exists! Storing old file as safeguard at: {newDest}")  # Prints out a heads up to user that backup has been created

    # Opens (or creates) file for storing data as a binary number
    dbfile = open(dest, 'wb') 

    # Dumps data in file
    pickle.dump(data, dbfile)

    # Closes file
    dbfile.close()



#################################
# --- Data loader function --- #
#################################

def loadData(target : str, showKeys = False):
    # Opens target file in bindary reading mode
    dbfile = open(target, 'rb')

    # Reads data
    db = pickle.load(dbfile)

    if showKeys:
        for key in db:
            print(key, '=>', db[key])
    
    # Closes file
    dbfile.close()

    return db



########################
# --- Test program --- #
########################

if __name__ == "__main__":
    cwd = os.getcwd() 
    print(cwd)

    data1 = [1, 2, 3]
    data2 = ['a', 'b', 'c']

    data = {"data1" : data1, "data2" : data2}

    storeData(data, cwd, "db")
    loadData("db", showKeys = True)
