import ElectroPredictor as ep
import time

print('Welcome to ElectroPredictor\n')

filename = ep.getSDFfile()

start_time = time.time()
print(f'Starting ElectroPredictor for {filename}...\n')

ep.createDirectories()
ep.createFiles(filename)
ep.predictElectrophilicity(filename)
ep.deleteFiles()

finishTime = round(((time.time() - start_time)/60.0),3)

print("\n\nElectrophiliciy calculations took %s minutes" % finishTime)

ep.jvm.stop()

